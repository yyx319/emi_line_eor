# To start a new track, this region should be set to a very specific a constraining
# volume containing EXCLUSIVELY the galaxy of interest.
module load intel-fortran/16.0.3
module load intel/2019u5
module load openmpi/intel/11.1/1.4.2


# input
# 1: simulation name 
# 2: initialize or track
echo $1
echo $2
track_dir=/home/yy503/Desktop/simulation_code/ramses_tools/tracker
sim_dir=/data/ERCblackholes2/smartin/Dwarf1
pp_dir=/data/ERCblackholes3/yuxuan/dwarf_data/post_processing

if [[ $2 == 'init' ]]; then
    xmi=0.50
    xma=0.52 
    ymi=0.47 
    yma=0.49
    zmi=0.48 
    zma=0.50
    reg="-xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma"

    op_idx=30

    
    output_track="$pp_dir/$1/track_000$op_idx.txt"

    inout="-inp $sim_dir/$1/output_000$op_idx -otk $output_track" 
    opts="$reg $inout"
    echo "------------------------"
    echo "particle_tracker.out $opts"
    echo "------------------------"
    echo 
    $track_dir/particle_tracker.out $opts
fi



if [[ $2 == 'track' ]]; then

for op_idx in {51..60..1}
  do
    op_idx_pre=$(($op_idx-1))
    echo $op_idx
    # Now we can switch to a more global tracking within a larger volume. Disabling tracking in
    # the entire simulation is a better practice, as it makes the code significantly
    # faster (it skips particles), but it still works equally well, as the only requirement
    # is finding the actual particles
    mival=0.45
    maval=0.55
    reg="-xmi $mival -xma $maval -ymi $mival -yma $maval -zmi $mival -zma $maval"

    # For the next iteration, we now use our previous tracked particles around the previous center
    # to find the new position of the galaxy
    input_track="$pp_dir/$1/track_000$op_idx_pre.txt"
    output_track="$pp_dir/$1/track_000$op_idx.txt"
    inout="-inp $sim_dir/$1/output_000$op_idx -new 1 -trk 1 -trf $input_track -otk $output_track"
    opts="$reg $inout"
    echo "------------------------"
    echo "particle_tracker.out $opts"
    echo "------------------------"
    echo 
    $track_dir/particle_tracker.out $opts

  done
fi