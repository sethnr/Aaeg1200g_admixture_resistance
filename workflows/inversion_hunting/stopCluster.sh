#!/bin/bash

export CLUSTER_NAME="redmond-aaa1-cluster1";
export ZONE="australia-southeast2";

gcloud container clusters delete $CLUSTER_NAME  --zone $ZONE
