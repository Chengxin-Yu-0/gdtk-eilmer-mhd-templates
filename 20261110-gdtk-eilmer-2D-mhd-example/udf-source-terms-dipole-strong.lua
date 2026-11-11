-- Source Terms Calculation Function - By Daniel Lewis, Alexis Lefevre, Sebastiaan van Oeveren, Aaron Cunningham and Peter Jacobs

-- Edits by Chengxin Yu, 11/11/2026
-- Cleaned up some outdated lines
-- Added more comments
-- Note that Raiser's conductivity formulae for ONLY strongly ionised plasma is included
-- This is valid as long as mole fraction of e- is greater than 1e-3
-- This copy is specifically for a magnetic field modelled by dipole equations

-----------------------------------------------------------------------------------------------------------------------------------------------
--CONSTANTS
-----------------------------------------------------------------------------------------------------------------------------------------------

-- some constants needed for the calculations
local mass_e     = 9.10938356e-31 -- mass of electron [kg]
local Na         = 6.02e23        -- Avogadro's constant [mol^-1]
local n_e_min    = 1.0e19         -- minimum electron number density for Lorentz force calculation [m^-3], this was set by Seb to avoid numerical issues in low density/freestream flow 
local Mr_e       = 5.48579903e-07 -- molecular mass of electron [kg/mol]

-----------------------------------------------------------------------------------------------------------------------------------------------
--MAGNEITC FIELD (DIPOLE EQUATION)
-----------------------------------------------------------------------------------------------------------------------------------------------

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

-----------------------------------------------------------------------------------------------------------------------------------------------
--CONDUCTIVITY FUNCTION
-----------------------------------------------------------------------------------------------------------------------------------------------

function conductivity(n_e, T)
    local sigma = 0.0

    if n_e < n_e_min then
        -- if the electron number density is below the minimum threshold, set conductivity to zero
        sigma    = 0.0 
    else
        -- setting up variables for Raiser's conductivity calculation
        local Te        = 0.0
        local lnA       = 0.0
        local K_to_eV   = 11604.525  -- conversion between Kelvin and eV [K eV-1]

        Te        = T/K_to_eV
        lnA       = 13.57+1.5*math.log10(Te)-0.5*math.log10((n_e+1.0)*1.0e-6)
        sigma     = (1.9e2*(Te^(1.5))/lnA)*100.0  -- [Ohm-1 m-1]
    end

    -- This is just to handle some extreme cases that can occur in low density/freestream flow
    if sigma == nil or sigma > 1.0e5 then 
        sigma = 0.0
    end

    return sigma
end -- conductivity

-----------------------------------------------------------------------------------------------------------------------------------------------
--LORENTZ SOURCE FUNCTION
-----------------------------------------------------------------------------------------------------------------------------------------------

function Lorentz_src(Bx,By,Bz,ux,uy,uz,sigma)
    local F = {}
    F.x = sigma*((uz*Bx-ux*Bz)*Bz-(ux*By-uy*Bx)*By)
    F.y = sigma*((ux*By-uy*Bx)*Bx-(uy*Bz-uz*By)*Bz)
    --F.z = sigma*((uy*Bz-vz*By)*By-(uz*Bx-ux*Bz)*Bx)  --uncomment for 3D
    F.z = 0.0
    return F
end --Lorentz_src

-----------------------------------------------------------------------------------------------------------------------------------------------
--SOURCE TERMS FUNCTION
-- computes the Lorentz force source term which appears in the resistive MHD equations
-----------------------------------------------------------------------------------------------------------------------------------------------

function sourceTerms(t, cell)
   
    -- import cell data to local variables
    src = {}
    local ux            = cell.velx
    local uy            = cell.vely  
    local uz            = cell.velz  
    local massf         = cell.massf
    local rho           = cell.rho
    local T             = cell.T

     -- calculating local electron number density, n_e, from mass fraction of e-, density, and relative molecular mass of e-
    local n_e = 0.0
    if cell.massf["e-"] == nil then 
        n_e = 0.0
    else 
        n_e = Na*cell.massf["e-"]*cell.rho/Mr_e
    end

    if n_e < n_e_min then
        src.momentum_x   = 0.0
        src.momentum_y   = 0.0
        src.momentum_z   = 0.0
        src.total_energy = 0.0
    else
         -- this is assuming magnetic field must be calculated using dipole equations from its cell-centred coordinate at each cell
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
