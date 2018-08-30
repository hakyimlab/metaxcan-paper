#!/usr/bin/env bash

#save data to s3
tar --dereference -czvpf data.tar.gz data/
aws s3 cp data.tar.gz s3://imlab-open/Data/MetaXcan/paper_data/metaxcan_paper_support_data_rev4.tar.gz
rm data.tar.gz