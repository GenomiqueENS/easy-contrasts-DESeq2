#!/bin/bash

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
SCRIPT_DIRNAME=$(basename "$SCRIPT_DIR")
DOCKER_IMAGE=${1:-genomicpariscentre/easycontrasts:2.0}
PROJECT_NAME=GSE107401
DATA_DIR=project_GSE107401
CORRESP_PATH=$DATA_DIR/ensembl_to_symbols.tsv

#
# Test methods
#

exit_if_fail() {

    if [ "$?" -ne 0 ]; then
        echo "ERROR: $1" >&2
        exit 1
    fi
}

is_file_exists() {

    if [ ! -s "$1" ]; then
        echo "ERROR: expected $1 file is missing." >&2
        exit 1
    fi
}

test_deseq2_comparison_files() {
    f="${PREFIX}_${PROJECT_NAME}-diffana_$1.tsv"
    is_file_exists "$f"
    logfc_col_name=$(head -n 1 "$f" | cut -f 4 | sed 's/log2foldchange //')
    expected_col_name=$(echo $1 | tr '_' ' ')
    if [ "$logfc_col_name" != "$expected_col_name" ]; then
        echo "ERROR: unexpected column name \"log2foldchange $logfc_col_name\" in $f file." >&2
        exit 1
    fi
}

test_deseq2_common_output() {

    is_file_exists "${PREFIX}_${PROJECT_NAME}-processed_dds_object.rds"

    # Test if normalisation files has been generated
    for s in $(echo normalisation_rawCountMatrix normalisation_normalisedCountMatrix); do
        for e in $(echo tsv rds); do
            is_file_exists "${PREFIX}_${PROJECT_NAME}-$s.$e"
        done
    done
}

test_deseq2_complex_output() {

    # Test if comparison output file has been generated
    for c in $(cut -f 1 ../project_GSE107401/deseq2_GSE107401-comparisonFile.txt); do
        test_deseq2_comparison_files "$c"
    done
}

#
# Tests
#

complex_mode_test() {
    echo "* Test complex mode"

    DESEQ_MODEL="~Condition+FooBar+Condition:FooBar"
    DESIGN_PATH=$DATA_DIR/deseq2_GSE107401-deseq2Design-complex.txt
    COMPARISON_PATH=$DATA_DIR/deseq2_GSE107401-comparisonFile.txt
    PREFIX=deseq2-complex
    RESULT_LINE_COUNT=2
    declare -A RESULT_LINE_COUNT_DICT
    RESULT_LINE_COUNT_DICT["KO_vs_WT"]=90
    RESULT_LINE_COUNT_DICT["Foo_vs_Bar"]=2
    RESULT_LINE_COUNT_DICT["FooKO_vs_FooWT"]=12

    # Remove previous output files
    rm -rf ${PREFIX}_* 2> /dev/null

    docker run \
    --rm \
    -v $(readlink -f ..):$(readlink -f ..) \
    -w $(readlink -f ..) \
    -u $(id -u):$(id -g) \
    $DOCKER_IMAGE \
    Rscript -e "rmarkdown::render(
        input = '01_normDiffana.Rmd',
        output_file = './$SCRIPT_DIRNAME/${PREFIX}_$PROJECT_NAME.html',
        params = list(projectName = '$PROJECT_NAME',
                      designPath = './$DESIGN_PATH',
                      comparisonPath = './$COMPARISON_PATH',
                      correspPath = './$CORRESP_PATH',
                      deseqModel = '$DESEQ_MODEL',
                      prefix = './$SCRIPT_DIRNAME/${PREFIX}_'))" > /dev/null
    exit_if_fail "Fail to execute R script with Docker."

    # Test output files
    test_deseq2_common_output
    test_deseq2_complex_output

    # Compare DESeq2 result files
    for c in $(cut -f 1 ../project_GSE107401/deseq2_GSE107401-comparisonFile.txt); do
        ./compare-deseq2-output.py --line-count $RESULT_LINE_COUNT ${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv expected-v1/${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv
        exit_if_fail "DESeq2 output files comparison failed."
        ./compare-deseq2-output.py --line-count $RESULT_LINE_COUNT ${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv expected-v2/${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv
        exit_if_fail "DESeq2 output files comparison failed."
        if [[ -v RESULT_LINE_COUNT_DICT[$c] ]]; then
            ./compare-deseq2-output.py --line-count ${RESULT_LINE_COUNT_DICT[$c]} ${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv expected-v1/${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv
            exit_if_fail "DESeq2 output files comparison failed."
            ./compare-deseq2-output.py --line-count ${RESULT_LINE_COUNT_DICT[$c]} ${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv expected-v2/${PREFIX}_${PROJECT_NAME}-diffana_$c.tsv
            exit_if_fail "DESeq2 output files comparison failed."
        fi
    done

    # Test if HTML output file exists
    is_file_exists ${PREFIX}_${PROJECT_NAME}.html
}

one_reference_mode_test() {
    echo "* Test one reference mode"

    DESEQ_MODEL="~Condition"
    DESIGN_PATH=$DATA_DIR/deseq2_GSE107401-deseq2Design-reference.txt
    COMPARISON_PATH=
    PREFIX=deseq2-reference

    # Remove previous output files
    rm -rf ${PREFIX}_* 2> /dev/null

    docker run \
    --rm \
    -v $(readlink -f ..):$(readlink -f ..) \
    -w $(readlink -f ..) \
    -u $(id -u):$(id -g) \
    $DOCKER_IMAGE \
    Rscript -e "rmarkdown::render(
        input = '01_normDiffana.Rmd',
        output_file = './$SCRIPT_DIRNAME/${PREFIX}_$PROJECT_NAME.html',
        params = list(projectName = '$PROJECT_NAME',
                    designPath = './$DESIGN_PATH',
                    correspPath = './$CORRESP_PATH',
                    deseqModel = '$DESEQ_MODEL',
                    prefix = './$SCRIPT_DIRNAME/${PREFIX}_'))" > /dev/null
    exit_if_fail "Fail to execute R script with Docker."

    # Test output files
    COMPARISON="KO_vs_WT"
    test_deseq2_common_output
    test_deseq2_comparison_files $COMPARISON

    # Compare DESeq2 result files
    for i in $(echo 1 10 100); do
        ./compare-deseq2-output.py --line-count $i ${PREFIX}_${PROJECT_NAME}-diffana_$COMPARISON.tsv expected-v1/${PREFIX}_${PROJECT_NAME}-diffana_$COMPARISON.tsv
        exit_if_fail "DESeq2 output files comparison failed."
        ./compare-deseq2-output.py --line-count $i ${PREFIX}_${PROJECT_NAME}-diffana_$COMPARISON.tsv expected-v2/${PREFIX}_${PROJECT_NAME}-diffana_$COMPARISON.tsv
        exit_if_fail "DESeq2 output files comparison failed."
    done

    # Test if HTML output file exists
    is_file_exists ${PREFIX}_${PROJECT_NAME}.html
}

multiple_references_mode_test() {
    echo "* Test one multiple references mode"

    DESEQ_MODEL="~Condition"
    DESIGN_PATH=$DATA_DIR/deseq2_GSE107401-deseq2Design-multi-references.txt
    COMPARISON_PATH=
    PREFIX=deseq2-multi-references
    declare -A RESULT_LINE_COUNT_DICT
    RESULT_LINE_COUNT_DICT["R0_vs_R1"]=2
    RESULT_LINE_COUNT_DICT["R0_vs_R2"]=3
    RESULT_LINE_COUNT_DICT["R2_vs_R1"]=1

    # Remove previous output files
    rm -rf ${PREFIX}_* 2> /dev/null

    docker run \
    --rm \
    -v $(readlink -f ..):$(readlink -f ..) \
    -w $(readlink -f ..) \
    -u $(id -u):$(id -g) \
    $DOCKER_IMAGE \
    Rscript -e "rmarkdown::render(
        input = '01_normDiffana.Rmd',
        output_file = './$SCRIPT_DIRNAME/${PREFIX}_$PROJECT_NAME.html',
        params = list(projectName = '$PROJECT_NAME',
                      designPath = './$DESIGN_PATH',
                      correspPath = './$CORRESP_PATH',
                      deseqModel = '$DESEQ_MODEL',
                      prefix = './$SCRIPT_DIRNAME/${PREFIX}_'))" > /dev/null
    exit_if_fail "Fail to execute R script with Docker"

    # Test output files
    test_deseq2_common_output

    for c in $(echo "R0_vs_R1" "R0_vs_R2" "R2_vs_R1"); do
        test_deseq2_comparison_files $c

        # Compare DESeq2 result files
        ./compare-deseq2-output.py --line-count ${RESULT_LINE_COUNT_DICT[$c]} ${PREFIX}_${PROJECT_NAME}-diffana_${c}.tsv expected-v1/${PREFIX}_${PROJECT_NAME}-diffana_${c}.tsv
        exit_if_fail "DESeq2 output files comparison failed."
        #./compare-deseq2-output.py --line-count ${RESULT_LINE_COUNT_DICT[$c]} ${PREFIX}_${PROJECT_NAME}-diffana_${c}.tsv expected-v2/${PREFIX}_${PROJECT_NAME}-diffana_$COMPARISON.tsv
        #exit_if_fail "DESeq2 output files comparison failed."
    done

    # Test if HTML output file exists
    is_file_exists ${PREFIX}_${PROJECT_NAME}.html
}

#
# Main
#

# Go to the script directory
cd "$SCRIPT_DIR"

complex_mode_test
one_reference_mode_test
multiple_references_mode_test

exit 0
