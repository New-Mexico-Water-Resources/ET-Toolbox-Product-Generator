from ETtoolbox.ETtoolbox_coarse import ET_toolbox_hindcast_tile

working_directory = "~/data/coarse_NRT_testing"
static_directory = "~/data/PTJPL_static"
tile = "11SPS"

ET_toolbox_hindcast_tile(tile=tile, working_directory=working_directory, static_directory=static_directory)