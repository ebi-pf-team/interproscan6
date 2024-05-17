# Linting

The `nf-core` suite includes its own linter. However, only `nf-core` packages can use the `nf-core` linter. While `InterProScan6` is not integrated into `nf-core` we use the `linter-rules-for-nextflow` linter from [AWSlabs](https://github.com/awslabs/linter-rules-for-nextflow).

## Setup

You can build the linter mannually by following the instructions [here](https://github.com/awslabs/linter-rules-for-nextflow?tab=readme-ov-file#build), alternatively (and recommend) to build the docker image for `linter-rules-for-nextflow` following [these instructions](https://github.com/awslabs/linter-rules-for-nextflow?tab=readme-ov-file#docker).

_Run all commands from the root of the repo._

1. Clone the `linter-rules-for-nextflow-repo`

```bash
git clone https://github.com/awslabs/linter-rules-for-nextflow.git
```

2. Use the bash script `drop_rules.sh` to drop rules we don't want/need in `InterProScan6`: 

* Add the corresponding XML line from the `general.xml` ruleset file from `linter-rules-for-nextflow` to the `RULES_TO_REMOVE` variable in `drop_rules.sh`:

```bash
RULES_TO_REMOVE=(
    '<rule class="software.amazon.nextflow.rules.ModuleIncludedTwiceRule"/>'
    '<rule class="rule.to.remove.dont.forget.to.escape.the.termianl.>"/>'
)
```
> https://github.com/awslabs/linter-rules-for-nextflow/blob/main/linter-rules/src/main/resources/rulesets/general.xml

* Add the files for the rule to the `FILES_TO_REMOVE` array in `drop_rules.sh`

```bash
# $1 is the path to the linter-for-nextflow dir
FILES_TO_REMOVE=(
    "$1/linter-rules/src/main/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRule.groovy"
    "$1/linter-rules/src/test/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRuleTest.groovy"
)
```

_To find the file paths, search the GitHub `linter-for-nextflow` repo with the name of the rule, e.g. "ModuleIncludedTwiceRule"._

* Run `drop_rules.sh`

```bash
# from the root of the IPS6 repo
bash linter/drop_rules.sh \
    <path to local copy of linter-rules-for-nextflow>
```

* Build the docker image

```bash
cd linter-rules-for-nextflow
docker build -t -linter-rules-for-nextflow .
cd ..  # return too interproscan6 root
```

## Run linting

Point the terminal at the root of the `InterProScan6` repository, and run this command for linting with the general rules:

```
docker run -v $PWD:/data -e ruleset=general linter-rules-for-nextflow
```

The container is configured to (by default) run all rules herein against all Nextflow (`*.nf`) files found in the data volume. 
