#!/bin/bash

# drop_rules.sh
# Drop linting rules we don't need/want for InterProScan6
# $1 = path to local copy of linter-rules-for-nextflow

RULES_TO_REMOVE=(
    '<rule class="software.amazon.nextflow.rules.ModuleIncludedTwiceRule"/>'
)
FILES_TO_REMOVE=(
    "$1/linter-rules/src/main/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRule.groovy"
    "$1/linter-rules/src/test/groovy/software/amazon/nextflow/rules/ModuleIncludedTwiceRuleTest.groovy"
)

RULES_XML_PATH="$1/linter-rules/src/main/resources/rulesets/general.xml"

# Check if RULES_XML_PATH exists
if [ ! -f "$RULES_XML_PATH" ]; then
    echo "Error: linter-for-nextflow XML file does not exist at $RULES_XML_PATH"
    exit 1
fi

for RULE in "${RULES_TO_REMOVE[@]}"
do
    echo "Removing the rule: $RULE"
    ESCAPED_STRING=$(sed 's/[&/\]/\\&/g' <<< "$RULE")
    TEMP_FILE=$(mktemp)
    sed "/$ESCAPED_STRING/d" "$RULES_XML_PATH" > "$TEMP_FILE"
    mv "$TEMP_FILE" "$RULES_XML_PATH"
done

for FILE in "${FILES_TO_REMOVE[@]}"
do
    echo "Removing file $FILE"
    rm -rf "$FILE"
done
