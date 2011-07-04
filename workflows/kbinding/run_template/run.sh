#!/bin/sh

BASE_DIR=`pwd`

. /etc/bashrc
module purge
module load intel openmpi
module load gcc/4.4.0
module load python
module load vmd
source /home/dacaplan/ENV/bin/activate
MPIFLAGS="-mca btl_sm_num_fifos 7 -mca mpi_paffinity_alone 1 -np $(wc -l $PBS_NODEFILE | gawk '{print $1}')"

SHM_DIR=/dev/shm/dacaplan
MPIRUN=/scinet/gpc/mpi/openmpi/1.4.1-intel-v11.0-ofed/bin/mpirun
NAMD=~/bin/namd2

STARTING_PDB=$BASE_DIR/md.pdb
PDB=$BASE_DIR/$coordinates
PSF=$BASE_DIR/md.psf

OUTPUT_DIR=$BASE_DIR/$id
SCRATCH_DIR=`echo $BASE_DIR | sed -e 's/\/project\/pomes\/dacaplan\//\/scratch\/dacaplan\/scratch\//g'`
SCRATCH_DIR=$SCRATCH_DIR/$id

# list of dirs to create if they don't exist
DIRS="$OUTPUT_DIR $SHM_DIR $SCRATCH_DIR"
# files to copy to working directory
COPY_FILES="$PSF $BASE_DIR/*.tcl $BASE_DIR/extrabonds $OUTPUT_DIR/restart.tar.gz"
    
# Create required directories
for d in $DIRS; do 
    if [ ! -e "$d" ]; then
        echo "Creating dir $d"
        mkdir -p $d
    fi
done

for f in $COPY_FILES; do
    cp $f $SHM_DIR/
done

echo "$sequence" > $OUTPUT_DIR/CURRENT
# set up the namd input file
cat $BASE_DIR/md.setup.conf | sed -e "s/__K__/$force/" | sed -e "s/__COORD__/$coordinate/" > $SHM_DIR/md.setup.conf
cat $BASE_DIR/md.restart.conf | sed -e "s/__K__/$force/" | sed -e "s/__COORD__/$coordinate/" > $SHM_DIR/md.restart.conf

cd $SHM_DIR

# IF restart files don't exist then we have to start from scratch
if [ ! -e "restart.tar.gz" ]; then
    echo "Restart files not found for replica $id"
    # if the PDB file doesn't exist then stage it
    if [ ! -e "${PDB}.gz" ]; then
        echo "PDB file ${PDB}.gz not found"
        # stage it
        python /home/dacaplan/projects/labwork/workflows/kbinding/replica_creation/stage.py $STARTING_PDB $template $PDB
        gzip $PDB
        if [ ! -e "${PDB}.gz" ]; then
            echo "Error staging the PDB file ${PDB}.gz from template $template and starting pdb $STARTING_PDB"
            exit 1
        fi
    fi
    cp ${PDB}.gz $SHM_DIR/md.pdb.gz
    gunzip $SHM_DIR/md.pdb.gz
    # now run md.setup.conf (short mini+equil to generate restart files)
    $MPIRUN $MPIFLAGS $NAMD $SHM_DIR/md.setup.conf >> $SHM_DIR/md.setup.log
    cp $SHM_DIR/md.setup.log $OUTPUT_DIR/
    if [ ! -e "restart.coor" ]; then
        echo "Error with initial MD step, check $OUTPUT_DIR/md.setup.log"
        exit 1
    fi
else
    # untar restart
    tar xzvf $SHM_DIR/restart.tar.gz
    rm $SHM_DIR/restart.tar.gz
    # copy starting pdb file 
    cp ${PDB}.gz $SHM_DIR/md.pdb.gz
    gunzip $SHM_DIR/md.pdb.gz
fi

RETURN_CODE=0

function finish {
    printf "Saving results from directory $SHM_DIR to $OUTPUT_DIR on "; date

    CORRUPT=0
    if [ -f "restart.vel" ] && [ -f "restart.vel.old" ]; then
        S=$(stat -c%s "restart.vel")
        S_OLD=$(stat -c%s "restart.vel.old")
        if [ "$S" -ne "$S_OLD" ]; then
            CORRUPT=1
        fi
    fi

    if [ -f "restart.coor" ] && [ -f "restart.coor.old" ]; then
        S=$(stat -c%s "restart.coor")
        S_OLD=$(stat -c%s "restart.coor.old")
        if [ "$S" -ne "$S_OLD" ]; then
            CORRUPT=1
        fi
    fi

    if [ "$CORRUPT" -eq "1" ]; then
        echo "restart files corrupt, replacing with old versions..."
        mv restart.xsc.old restart.xsc
        mv restart.coor.old restart.coor
        mv restart.vel.old restart.vel
    fi

    if [ ! -e "md.dcd" ]; then
        echo "No DCD file found... something is wrong!"
        gzip md.log
        mv md.log.gz $OUTPUT_DIR/md.error.log.gz
        RETURN_CODE=1
    else
        rm restart.*.old
        rm $SHM_DIR/*.BAK
        # if we are sure that things ran properly
        gzip md.log
        mv md.log.gz $OUTPUT_DIR/
        tar zcf $OUTPUT_DIR/restart.tar.gz restart.*

        FRAMEINDEX=`~/bin/catdcd md.dcd | grep Total | cut -d ":" -f 2 | sed -e 's/ //g'`
        ~/bin/catdcd -o $sequence.pdb -otype pdb -stype pdb -s md.pdb -first $FRAMEINDEX -last $FRAMEINDEX md.dcd 2>&1 >> /dev/null
        gzip $sequence.pdb
        mv $sequence.pdb.gz $SCRATCH_DIR/
        #cp $SHM_DIR/md.dcd $SCRATCH_DIR/$sequence.dcd
        ANALYSIS=/home/dacaplan/projects/labwork/workflows/kbinding/analysis.py
        python $ANALYSIS -p md.psf md.dcd >> $OUTPUT_DIR/distances
    fi

    printf "Cleaning up temporary directory $SHM_DIR at "; date
    rm -rf $SHM_DIR
    printf "done at "; date
}

function trap_term {
    printf "Trapped term (soft kill) signal on "; date
    gzip md.log
    mv md.log.gz $OUTPUT_DIR/md.error.log.gz
    exit 1
}

trap "trap_term" TERM

printf "Starting run at "; date
echo "Running: $MPIRUN $MPIFLAGS $NAMD $SHM_DIR/md.restart.conf >> $SHM_DIR/md.log"
if [ "$MPIRUN" == "" ]; then
    $NAMD +p4 $SHM_DIR/md.restart.conf >> $SHM_DIR/md.log
else
    $MPIRUN $MPIFLAGS $NAMD $SHM_DIR/md.restart.conf >> $SHM_DIR/md.log
fi
printf "Run ended cleanly at "; date
finish
exit $RETURN_CODE
