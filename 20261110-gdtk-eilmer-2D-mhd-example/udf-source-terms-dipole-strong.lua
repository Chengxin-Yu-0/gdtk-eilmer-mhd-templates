-- udf-source-terms.lua
-- Source term for including Lorentz Force as a body force in the flow
-- original by Alexis Lefevre - May 2019
-- revisions by Daniel Smith 7-May-2020
--      changed magnetic field formulation to be correct for spherical magnet, rather than using theoretical dipole, equation from http://farside.ph.utexas.edu/teaching/jk1/lectures/node61.html
-- Note this script does not account for Hall effects, nor does it account for the MHD source term in the energy conservation equation. 
-- Electrical conductivity is a scalar value
-- major revisions by Aaron Cunningham May 2023
--      updated to pull in cell data for magnetic fields from flow data.  mag field inputs to flow field is from udf-supervisor-file
--      sigmoid scale added instead of schedule table.
--      updated for 3D B field
-- Cleaning up by Chengxin
--      removed magnetic field ramping

--some convenient functions for use later
local exp   = math.exp
local sqrt  = math.sqrt
local pi    = math.pi
local log   = math.log

-- some constants needed for the calculations
local W_e         = 5.485799e-7    -- Units are kg/mol
local mass_e      = 9.10938356e-31 -- mass of electron in kg
local Na          = 6.02e23        --Avogadro's constant
local n_e_min     = 1.0e19

function magnetic_field(t, x, y)
  local B_vector = {}

  local B0 = 0.8
  local Rb = 0.01905
  local xc = x - Rb
  local yc = y 
  local r = sqrt(xc^2 + yc^2)

  local Bx = -((B0*Rb^3)/(2*r^5)) * (2*(xc)^2 - (yc)^2) -- Tesla
  local By = -((B0*Rb^3)/(2*r^5)) * (3*(xc)*(yc))  -- Tesla

  B_vector.x = Bx
  B_vector.y = By
  B_vector.z = 0.0

  return B_vector
end

function conductivity(n_e, T)
    local sigma = 0.0

    --electrical conductivity from Raiser for a strongly ionized plasma
    if n_e < n_e_min then
        sigma    = 0.0
    else
        local Te  = 0.0; 
        local lnA = 0.0
        local K_to_eV   = 11604.525  --conversion between Kelvin and eV  --K eV-1
        
        Te        = T / K_to_eV    --acc T is defined in sourceTerms - added to function args
        lnA       = 13.57 + 1.5 * math.log10(Te) - 0.5 * math.log10((n_e+1.0)*1.0e-6)
        sigma     = (1.9e2 * (Te^(1.5)) / lnA) * 100.0  --Ohm-1 m-1
    end

    if sigma == nil or sigma > 1.0e5 then  -- This is just to handle some extreme cases that can occur in low density/freestream flow
        sigma = 0.0
    end
    return sigma -- conductivity
end

function Lorentz_src(Bx,By,Bz,ux,uy,uz,sigma)
    local F = {}
    F.x = sigma*((uz*Bx-ux*Bz)*Bz-(ux*By-uy*Bx)*By)
    F.y = sigma*((ux*By-uy*Bx)*Bx-(uy*Bz-uz*By)*Bz)
    --F.z = sigma*((uy*Bz-vz*By)*By-(uz*Bx-ux*Bz)*Bx)  --uncomment for 3D
    F.z = 0.0
    return F
end --Lorentz_src

function sourceTerms(t, cell)
   
    -- computes the Lorentz force source term which appears in the resistive MHD equations

    -- import cell data to local variables
    local src = {}
    local ux            = cell.velx
    local uy            = cell.vely  
    local uz            = cell.velz

    local massf_e       = cell.massf["e-"]    
    local rho           = cell.rho
    local T             = cell.T
    local p             = cell.p

    --electron mass fraction
    local n_e = 0.0
    if cell.massf["e-"] == nil then 
        n_e = 0.0
    else 
        n_e = Na*cell.massf["e-"]*rho/W_e      --electron number density
    end

    if n_e < n_e_min then

        src.momentum_x   = 0.0 --   -sigma*By*(ux*By-uy*Bx)
        src.momentum_y   = 0.0 --    sigma*Bx*(ux*By-uy*Bx)
        src.momentum_z   = 0.0
        src.total_energy = 0.0
    
    else

        local x             = cell.x
        local y             = cell.y

        local B_vector = magnetic_field(t, x, y)
        local Bx = B_vector.x
        local By = B_vector.y
        local Bz = B_vector.z

        local sigma       = conductivity(n_e, T)                  
        local MHD_force   = Lorentz_src(Bx,By,Bz,ux,uy,uz,sigma)

        src.momentum_x   = MHD_force.x
        src.momentum_y   = MHD_force.y
        src.momentum_z   = 0.0
        src.total_energy = 0.0

    end
    return src
end --end function
