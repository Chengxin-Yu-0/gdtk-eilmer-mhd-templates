-- Source Terms Calculation Function - By Daniel Lewis, Alexis Lefevre, Sebastiaan van Oeveren, Aaron Cunningham and Peter Jacobs

-- Edits by Chengxin Yu, 28/10/2026
-- Cleaned up some outdated lines
-- Added more comments
-- Added Raizer's *interpolated* conductivity for 11 species air, previously only strongly ionised plasma was used
-- Note that Mr, mole fraction and absolute mole are hardcoded to 11 species air, this should be easy enough to change if needed for other gas models

-----------------------------------------------------------------------------------------------------------------------------------------------
--CONSTANTS
-----------------------------------------------------------------------------------------------------------------------------------------------

-- some constants needed for the calculations
local mass_e     = 9.10938356e-31 -- mass of electron [kg]
local Na         = 6.02e23        -- Avogadro's constant [mol^-1]
local n_e_min    = 1.0e19         -- minimum electron number density for Lorentz force calculation [m^-3], this was set by Seb to avoid numerical issues in low density/freestream flow 

-- relative molecular masses of atoms/molecules [kg/mol], referenced from air-11sp-gas-model.lua
-- for 11 species air: 'N2', 'O2', 'N', 'O', 'NO', 'N2+', 'O2+', 'N+', 'O+', 'NO+', 'e-'
local Mr_N2     = 2.80134000e-02
local Mr_O2     = 3.19988000e-02
local Mr_N      = 1.40067000e-02
local Mr_O      = 1.59994000e-02
local Mr_NO     = 3.00061000e-02
local Mr_N2p    = 2.80128514e-02
local Mr_O2p    = 3.19982514e-02
local Mr_Np     = 1.40061514e-02
local Mr_Op     = 1.59988514e-02
local Mr_NOp    = 3.00055514e-02
local Mr_e      = 5.48579903e-07

-----------------------------------------------------------------------------------------------------------------------------------------------
--CONDUCTIVITY FUNCTION
-----------------------------------------------------------------------------------------------------------------------------------------------

function conductivity(n_e,T,p,rho,vol,massf)
    sigma = 0.0
    
    -- total mass in a cell
    local mass = rho*vol

    -- calculating absolute mole of each species
    local mole_N2 = massf['N2']*mass/Mr_N2
    local mole_O2 = massf['O2']*mass/Mr_O2
    local mole_N = massf['N']*mass/Mr_N
    local mole_O = massf['O']*mass/Mr_O
    local mole_NO = massf['NO']*mass/Mr_NO
    local mole_N2p = massf['N2+']*mass/Mr_N2p
    local mole_O2p = massf['O2+']*mass/Mr_O2p
    local mole_Np = massf['N+']*mass/Mr_Np
    local mole_Op = massf['O+']*mass/Mr_Op
    local mole_NOp = massf['NO+']*mass/Mr_NOp
    local mole_e = massf['e-']*mass/Mr_e

    -- total absolute mole in a cell
    local mole = mole_N2+mole_O2+mole_N+mole_O+mole_NO+mole_N2p+mole_O2p+mole_Np+mole_Op+mole_NOp+mole_e

    -- mole fraction of e-
    local molef_e = mole_e/mole

    if n_e < n_e_min then
        -- if the electron number density is below the minimum threshold, set conductivity to zero
        sigma    = 0.0 
    else
        -- setting up variables for Raiser's conductivity calculation
        local Te        = 0.0
        local lnA       = 0.0
        local K_to_eV   = 11604.525  -- conversion between Kelvin and eV [K eV-1]

        -- electrical conductivity from Raiser for a strongly ionized plasma
        if molef_e >= 1e-3 then
            Te        = T/K_to_eV
            lnA       = 13.57+1.5*math.log10(Te)-0.5*math.log10((n_e+1.0)*1.0e-6)
            sigma     = (1.9e2*(Te^(1.5))/lnA)*100.0  -- [Ohm-1 m-1]

        -- electrical conductivity from Raiser for a weakly ionized plasma
        elseif molef_e <= 1e-4 then
            sigma = 2.82e-4*((n_e+1.0)*1.0e-6)/(5.3e9*p*133.322)*100.0  -- [Ohm-1 m-1]

        -- interpolating electrical conductivity
        else
            -- creating some temporary variables to hold the weak and strong conductivity values for interpolation
            local n_e_weak = 0.0
            local sigma_weak = 0.0
            local n_e_strong = 0.0
            local Te_temp = 0.0
            local lnA_temp = 0.0
            local sigma_strong = 0.0

            -- converting weak 1e-4 mole fraction of electron (1e-4*mole) to electron number density (n_e)
            n_e_weak        = (1e-4*mole)*Na/vol
            -- start of conductivity interpolation from weakly ionized plasma
            sigma_weak      = 2.82e-4*((n_e_weak+1.0)*1.0e-6)/(5.3e9*p*133.322)*100.0
            -- converting strong 1e-3 mole fraction of electron (1e-3*mole) to electron number density (n_e)
            n_e_strong      = (1e-3*mole)*Na/vol
            -- end of conductivity interpolation from strongly ionized plasma
            Te_temp         = T/K_to_eV
            lnA_temp        = 13.57+1.5*math.log10(Te_temp)-0.5*math.log10((n_e_strong+1.0)*1.0e-6)
            sigma_strong    = (1.9e2*(Te_temp^(1.5))/lnA_temp)*100.0 

            -- interpolate sigma at given molef['e-']
            local function interpolate(x1, y1, x2, y2, x)
                return y1+(x-x1)*(y2-y1)/(x2-x1)
            end

            sigma = interpolate(1e-4, sigma_weak, 1e-3, sigma_strong, molef_e)

        end
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
    F = {}
    F.x = sigma*((uz*Bx-ux*Bz)*Bz-(ux*By-uy*Bx)*By)
    F.y = sigma*((ux*By-uy*Bx)*Bx-(uy*Bz-uz*By)*Bz)
    --F.z = sigma*((uy*Bz-vz*By)*By-(uz*Bx-ux*Bz)*Bx)  -- uncomment for 3D
    F.z = 0
    return F
end -- Lorentz_src

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
    local vol           = cell.vol
    local p             = cell.p
    local T             = cell.T

    -- calculating local electron number density, n_e, from mass fraction of e-, density, and relative molecular mass of e-
    local n_e = 0.0
    if cell.massf["e-"] == nil then 
        n_e = 0.0
    else 
        n_e = Na*cell.massf["e-"]*cell.rho/Mr_e
    end

    -- Only calculate Lorentz force if sufficient ionization is reached, otherwise might get unstable/strange results
    if n_e < n_e_min then
        src.momentum_x   = 0.0
        src.momentum_y   = 0.0
        src.momentum_z   = 0.0
        src.total_energy = 0.0
    else

        -- this is assuming magnetic field has already been assigned to each cell, either through the udf supervisor file or some other method
        local Bx = cell.Bx
        local By = cell.By
        local Bz = cell.Bz

        sigma       = conductivity(n_e,T,p,rho,vol,massf)               
        MHD_force   = Lorentz_src(Bx,By,Bz,ux,uy,uz,sigma)

        src.momentum_x   = MHD_force.x
        src.momentum_y   = MHD_force.y
        src.momentum_z   = MHD_force.z
        src.total_energy = 0.0

    end
    return src
end --end function
