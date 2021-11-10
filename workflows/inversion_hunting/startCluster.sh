#!/bin/bash

export CLUSTER_NAME="redmond-aaa1-cluster1";
export ZONE="australia-southeast2";
#export MACHINE_TYPE="n1-standard-16"
export MACHINE_TYPE="n1-standard-8"

gcloud container clusters create $CLUSTER_NAME \
    --zone=$ZONE --num-nodes=1 \
    --machine-type=$MACHINE_TYPE \
    --scopes storage-rw \
    --image-type=UBUNTU \
    --disk-size=500GB \
    --enable-autoscaling \
    --max-nodes=96 \
    --min-nodes=0 \
    --service-account seth-redmond@slurm-cluster-308318.iam.gserviceaccount.com
