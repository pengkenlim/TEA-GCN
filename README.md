<img src="https://github.com/user-attachments/assets/682c7b18-fdc6-4353-bf35-28806b296484" alt="banner" width="400"/>

## TEA-GCN
### Two-Tier Ensemble Aggregation Gene Co-expression Network 

<img src="https://github.com/user-attachments/assets/f31ae18f-5846-49d7-b597-3234a7035ab2" alt="banner" width="700"/>

### What does this pipeline do?
This pipeline generates high-quality Gene Co-expression Networks (TEA-GCN ) that capture tissue/condition-specific co-expression.

### Step-by-step Guide
#### Step 1. Setting up
##### Clone repository to local machine
```
$ git clone https://github.com/pengkenlim/TEA-GCN.git
```
##### Create an environment, and install packages
```
$ cd TEA-GCN
$ virtualenv -p python3 .venv
$ source ./.venv/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
```

##### Activate environment for subsequent use
```
$ cd TEA-GCN
$ source ./.venv/bin/activate
```
