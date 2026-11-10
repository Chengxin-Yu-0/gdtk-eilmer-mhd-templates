-- Magnetic Field Input Function - By Aaron Cunningham 2023 with the help of Peter Jacobs
-- Function reads in the magnetic field_table.lua, and loads the magnetic fields to the flow data.
-- See field_table.lua included for example table format
-- Create_Grid_To_Ansys.py is example python script to interpolate field data from external source.

-- Edits by Chengxin Yu, 28/10/2026
-- Cleaned up some debugging lines

-----------------------------------------------------------------------------------------------------------------------------------------------
--CODE START
-----------------------------------------------------------------------------------------------------------------------------------------------

--load the field data from external source
dofile('field_table.lua')
num_blocks = how_many_blocks.num_blocks

print("\n LOADED SUPERVISOR \n")
print("blocks = ",num_blocks+1)

-----------------------------------------------------------------------------------------------------------------------------------------------
-- Entry point is atTimestepStart()
-- Note that atTimestepStart() is for transient solver
-----------------------------------------------------------------------------------------------------------------------------------------------
function atTimestepStart(tab)

    -- tab contains time, step, timeStep
    local step = tab.step

    -- only updates the B field after running the simulation, does not update upon prepping (snapshot 0000)
    -- if using not using scaler here set time steps to 3 or more to avoid calculation load.
    if step<5 then 
        
        -- user to add subfunctions if magnetic field needs to be incrementally ramped up to the desired field strength
        -- may be helpful for stronger magnetic fields
        local scaler = 1

        for block_index = 0,num_blocks,1 do
            blk_detail = infoFluidBlock(block_index)
            no_i = blk_detail.imax-blk_detail.imin
            no_j = blk_detail.jmax-blk_detail.jmin

            --add onebecause lua indexs arrays from 1
            num_cells = (no_i+1)*(no_j+1)

            -- march order - i1,j1 i2,j1 i3,j1...  i1,j2 i2,j2
            -- cell_index assumes above march order, ensure field_table is in the same order, should match the flow files
            cell_index = 0
            for cell_j = 1,no_j+1,1 do
                for cell_i = 1,no_i+1,1 do
                    cell_index = cell_index+1

                    -- in the example field_table.lua the Bx is at index 4 etc., adjust if your table is different
                    Bx = blocks[block_index][cell_index][4]*scaler
                    By = blocks[block_index][cell_index][5]*scaler
                    Bz = blocks[block_index][cell_index][6]*scaler

                    setBxyzInFluidCell(block_index, cell_i-1, cell_j-1, 0, Bx, By, Bz)     

                end --end cell_i loop
            end --end cell_j loop     
        end --end for 
    end --end if  
    return
end
