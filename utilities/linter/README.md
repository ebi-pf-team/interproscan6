# Linting

## Set up

Install using Git clone and the docker image from the `linter-rules-for-nextflow` repo

> https://github.com/awslabs/linter-rules-for-nextflow?tab=readme-ov-file#docker

## Drop rules we don't want/need

Use the bash script `drop_rules.sh` to drop rules we don't want/need in `InterProScan6`:

1. Add the corresponding XML line from the `general.xml` ruleset file from `linter-rules-for-nextflow` to the `RULES_TO_REMOVE` variable in `drop_rules.sh`:

```bash
RULES_TO_REMOVE=(
    '<rule class="software.amazon.nextflow.rules.ModuleIncludedTwiceRule"/>'
    '<rule class="rule.to.remove.dont.forget.to.escape.the.termianl.>"/>'
)
```

> https://github.com/awslabs/linter-rules-for-nextflow/blob/main/linter-rules/src/main/resources/rulesets/general.xml

3. Add the files for the rule to the `FILES_TO_REMOVE` array in `drop_rules.sh`

```bash
# $1 is the path to the linter-for-nextflow dir
FILES_TO_REMOVE=(
    "$1/linter-rules/src/main/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRule.groovy"
    "$1/linter-rules/src/test/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRuleTest.groovy"
)
```

_To find the file paths, search the GitHub `linter-for-nextflow` repo with the name of the rule, e.g. "ModuleIncludedTwiceRule"._

2. Run `drop_rules.sh`

```bash
# from the root of the IPS6 repo
bash linter/drop_rules.sh \
    <path to local copy of linter-rules-for-nextflow>
```

3. Rebuild the docker image

```bash
docker build -t -linter-rules-for-nextflow .
```
