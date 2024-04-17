# How to download ABIDE data in BIDS format?

## Create an AWS account
[Click here to find the tutorial](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html)

## Run AWS on your terminal
Before running aws s3 bucket commands the terminal, connect to your account. Open the AWS configuration file, stored locally on your machine. This file contains the sso_start_url, which you can copy & paste to a browser to connect to your AWS account.
```bash
#!/bin/bash
cat ~/.aws/config
```bash

### Once connected to your AWS Identity Account, open Access Keys. Copy and paste the AWS environmental variables to the bash terminal. Technically, you don't have to do it if you configure AWS by creating a credentials file.
 
## Run the AWS S3 bucket commands
### To list the datasets available on AWS s3
```bash
#!/bin/bash
aws s3 ls fcp-indi/data/Projects/ABIDE2/RawData/

### To download raw data from Georgetown University (GU_1) and Kennedy Krieger Institute (KKI_1)
```bash
#!/bin/bash
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-GU_1/ C:\LOCAL\PATH\GU_1
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-KKI_1/ C:\LOCAL\PATH\KKI_1
