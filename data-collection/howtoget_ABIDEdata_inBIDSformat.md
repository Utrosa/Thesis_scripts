# How to download ABIDE data in BIDS format?

## Create an AWS account
[Click here to read the setup guide](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html)

## Connect to AWS access portal
Before running AWS S3 bucket commands in the terminal, connect to your AWS access portal as an IAM user (not a root user). Open the AWS configuration file, stored locally on your machine. This file contains the sso_start_url, which you can copy & paste to a browser to connect to your AWS account.
```bash
#!/bin/bash
cat ~/.aws/config
```
Once connected to your account, open _Access Keys_. Then copy & paste the AWS environmental variables to the bash terminal. Technically, you don't have to do this, if you configure AWS by creating a credentials file. However, in this way, it will work for sure.
 
## Run the AWS S3 bucket commands
List the datasets available on AWS s3
```bash
#!/bin/bash
aws s3 ls fcp-indi/data/Projects/ABIDE2/RawData/
```
Download raw data from Georgetown University (GU_1) and Kennedy Krieger Institute (KKI_1)
```bash
#!/bin/bash
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-GU_1/ C:\LOCAL\PATH\GU_1
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-KKI_1/ C:\LOCAL\PATH\KKI_1
```