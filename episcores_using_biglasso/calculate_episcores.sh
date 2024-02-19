#!/bin/bash

prepData() {
    echo "Calculate EpiScores for a set of traits - data prep - procedure begins"
    SETTINGS=$1

    for TRAIT in $SETTINGS/*
    do
        echo "Trait: $TRAIT"
        
        Rscript --vanilla /Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/data_prep.R --settings $TRAIT
        # Rscript --vanilla ~/Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/data_prep.R --settings $TRAIT
		echo "Finished, check output size"
    done

    echo "Procedure finished."
}


trainEpiScores() {
    echo "Calculate EpiScores for a set of traits - models - procedure begins"
    SETTINGS=$1

    for TRAIT in $SETTINGS/*
    do
        echo "Trait: $TRAIT"
        
        Rscript --vanilla /Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/data_train.R --settings $TRAIT
		echo "Finished, check output size"
    done

    echo "Procedure finished."
}

testEpiScores() {
    echo "Calculate EpiScores for a set of traits - models - procedure begins"
    SETTINGS=$1

    for TRAIT in $SETTINGS/*
    do
        echo "Trait: $TRAIT"
        
        Rscript --vanilla /Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/data_test.R --settings $TRAIT
		echo "Finished, check output size"
    done

    echo "Procedure finished."
}

testEpiScores /Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/settings
