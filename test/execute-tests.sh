#!/bin/sh

DOCKER_IMAGE=genomicpariscentre/easycontrasts:2.0
PROJECT_NAME=GSE107401
DATA_DIR=project_GSE107401
CORRESP_PATH=$DATA_DIR/ensembl_to_symbols.tsv



#
# Test methods
#

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
# DESeq2 in complex mode
#

DESEQ_MODEL="~Condition+FooBar+Condition:FooBar"
DESIGN_PATH=$DATA_DIR/deseq2_GSE107401-deseq2Design-complex.txt
COMPARISON_PATH=$DATA_DIR/deseq2_GSE107401-comparisonFile.txt
PREFIX=deseq2-complex

# Remove existing output files
rm -rf ${PREFIX}_* 2> /dev/null

docker run \
-ti --rm \
-v $(readlink -f ..):$(readlink -f ..) \
-w $(readlink -f ..) \
-u $(id -u):$(id -g) \
$DOCKER_IMAGE \
Rscript -e "rmarkdown::render(
    input = '01_normDiffana.Rmd',
    output_file = '$PROJECT_NAME.html',
    params = list(projectName = '$PROJECT_NAME',
                  designPath = './$DESIGN_PATH',
                  comparisonPath = './$COMPARISON_PATH',
                  correspPath = './$CORRESP_PATH',
                  deseqModel = '$DESEQ_MODEL',
                  prefix = './test/${PREFIX}_'))"

# Test output files
test_deseq2_common_output
test_deseq2_complex_output

#
# DESeq2 in reference mode
#

DESEQ_MODEL="~Condition"
DESIGN_PATH=$DATA_DIR/deseq2_GSE107401-deseq2Design-reference.txt
COMPARISON_PATH=
PREFIX=deseq2-ref

# Remove output files
rm -rf ${PREFIX}_* 2> /dev/null


docker run \
-ti --rm \
-v $(readlink -f ..):$(readlink -f ..) \
-w $(readlink -f ..) \
-u $(id -u):$(id -g) \
$DOCKER_IMAGE \
Rscript -e "rmarkdown::render(
    input = '01_normDiffana.Rmd',
    output_file = '$PROJECT_NAME.html',
    params = list(projectName = '$PROJECT_NAME',
                  designPath = './$DESIGN_PATH',
                  correspPath = './$CORRESP_PATH',
                  deseqModel = '$DESEQ_MODEL',
                  prefix = './test/${PREFIX}_'))"

# Test output files
test_deseq2_common_output
test_deseq2_comparison_files "KO_vs_WT"

# Remove output files
rm -rf ${PREFIX}_* 2> /dev/null

exit 0
