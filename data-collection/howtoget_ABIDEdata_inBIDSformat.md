# How to download ABIDE data in BIDS format?

## Create an AWS account
[Click here to find the tutorial](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html)

## Run AWS on your terminal
### Before running aws s3 bucket commands the terminal, connect to your account. To do this, open the AWS configuration file, stored locally on your machine. This file contains the `sso_start_url`. Copy & paste it to a browser and connect to your AWS account with your username and password (not as as root user).
```
cat ~/.aws/config

### Once connected to your AWS Identity Account, open Access Keys. Select the tab corresponding to your system (Linux, PowerShell, ...). Copy and paste the AWS environmental variables to the terminal.
 
## Run the AWS S3 bucket commands
### To list the datasets available on AWS s3
```
aws s3 ls fcp-indi/data/Projects/ABIDE2/RawData/

### To download raw data from Georgetown University (GU_1), Kennedy Krieger Institute (KKI_1), and Oregon Health and Science University (OHSU) run the following commands in your terminal.
```
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-GU_1/ C:\LOCAL\PATH\GU_1
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-KKI_1/ C:\LOCAL\PATH\KKI_1
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ABIDEII-OHSU_1/ C:\LOCAL\PATH\OHSU_1