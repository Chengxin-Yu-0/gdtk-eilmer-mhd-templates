print("Sphere MHD with Dipole Magnetic Field")

--------------------------
----- Flow Condition -----
--------------------------

no_flow_times = 8.0

nsp, nmodes, gm = setGasModel('air-7sp-1T.gas')
config.reactions_file = 'air-7sp-9r-1T.chem'

--- X2 Zander Condition [Zander, 20??] ---
T_inf = 2182.866-- Kelvin  
p_inf= 870.9855 -- Pa, 1 atm
vel_inf = 8223.9  -- m/s
mf_inf = {['N2']=0.767, ['O2']=0.233}

-- this code block allows the magnetic field to be visualised upon prepping the sim
function initialState(x, y, z)
   B0 = 0.8
   Rb = 0.01905
   xc = x - Rb
   yc = y
   r = (xc^2 + yc^2)^0.5
   B0x = -((B0*Rb^3)/(2*r^5)) * (2*(xc)^2 - (yc)^2) -- Tesla
   B0y = -((B0*Rb^3)/(2*r^5)) * (3*(xc)*(yc))  -- Tesla
   -- visual check to ensure magnet is not within fluid grid
   if r < Rb then
      B0x = 600
      B0y = 600
   end
   return FlowState:new{p=p_inf, T=T_inf, velx=vel_inf, Bx=B0x, By=B0y, massf=massf_inf, T_modes={T_inf}}
end
initial = initialState
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf}, velx=vel_inf,massf=mf_inf}
            
flowDict = {
   initial=initial,
   inflow=inflow
}

bcDict = {
   inflow=InFlowBC_Supersonic:new{flowState=inflow},
   outflow=OutFlowBC_SimpleExtrapolate:new{},
   wall=WallBC_WithSlip0:new{}
}

makeFluidBlocks(bcDict, flowDict)
print('Fluid block constructed\n')

---------------------------
----- Config Settings -----
---------------------------

--source terms (for Lorentz force)
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms-dipole-strong.lua'
config.MHD = true

radius = 0.01905 -- consistent with 0.8T 38.1mm diameter sphere magnet

body_flow_time = radius/vel_inf
interpolation_delay = 4*body_flow_time
max_time = no_flow_times*body_flow_time
dt_plot = body_flow_time/5
print('==============================================')
print("body_flow_time      : ",body_flow_time)
print('Interpolation Delay : ',interpolation_delay)
print("Max Time            : ",max_time)
print("Snapshot interval   : ",dt_plot)
print('Total snapshots expected : ',max_time/dt_plot)
print('==================================')
config.solver_mode = 'transient'
config.reacting = true
config.flowstate_limits_max_temp = 50000.0
config.dimensions = 2
config.axisymmetric = true
config.viscous = true
config.flux_calculator = "hanel"
config.interpolation_order = 2
config.interpolation_delay = interpolation_delay
config.max_invalid_cells = 20
config.adjust_invalid_cell_data = true
config.max_time = max_time
config.max_step = 10000000
config.dt_init = 1.0e-10
config.cfl_value = 0.4
config.stringent_cfl = true
config.dt_plot = dt_plot
--