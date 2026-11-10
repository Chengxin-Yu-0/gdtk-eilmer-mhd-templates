
config.dimensions = 2

-----------------------
----- Create Grid -----
-----------------------

radius = 0.01905 -- consistent with 0.8T 38.1mm diameter sphere magnet
x_grid = 0.002 -- extend inflow grid by x_grid [mm]

A = Vector3:new{x=0.0, y=0.0} -- Stagnation point
B = Vector3:new{x=radius, y=radius} -- Shoulder
centre = Vector3:new{x=radius, y=0.0} -- Centre of the sphere/model

-- Some points to define the inflow and outflow domains
D_x = -x_grid
D_centre_x = 1.5*radius
C_x = B.x-1.5*x_grid

D = Vector3:new{x=D_x, y=0.0} -- Inflow edge
D_centre = Vector3:new{x=D_centre_x, y=0.0} -- Centre of the inflow arc
inflow_radius = math.abs(D_centre_x-D_x)
C_y = math.sqrt( radius^2 - (C_x-D_centre_x)^2 )
function circ_point(D_x, D_centre_x, C_x)
   radius = math.abs(D_centre_x-D_x)
   C_y = math.sqrt( radius^2 - (C_x-D_centre_x)^2 )
   return C_y
end
C_y = circ_point(D.x, D_centre.x, C_x)
C = Vector3:new{x=C_x, y=C_y} -- Outflow edge
print('Points Defined')

AB = Arc:new{p0=A, p1=B, centre=centre} -- Body
CB = Line:new{p0=C, p1=B} -- Outflow
-- DC = Bezier:new{points = {D, bez_D, bez_C, C}} -- Inflow
DC = Arc:new{p0=D, p1=C, centre=D_centre} -- Inflow
DA = Line:new{p0=D, p1=A} -- Symmetry
print('Lines Defined')

-- Number of cells
niv=40+1; njv=100+1;

quad1 = makePatch{north=CB, east=AB, south=DA, west=DC}
grid1 = StructuredGrid:new{psurface=quad1, niv=niv, njv=njv}
print('Grids Defined')

-- Fluidblocks
nib = 1; njb = 6;
print('Number of fluid blocks: ',nib*njb)
-- Define blocks
registerFluidGridArray{
   grid=grid1,
   nib=nib, njb=njb,
   fsTag="initial",
   shock_fitting=true,
   bcTags={west="inflow", east="wall", north="outflow"}
}

identifyGridConnections()
print('Grid Connection Identified')
