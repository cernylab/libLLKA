/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_minicif.h>

#include <stdio.h>
#include <stdlib.h>

static
void mirror(const LLKA_CifData *data)
{
    size_t blockIdx;

    for (blockIdx = 0; blockIdx < data->nBlocks; blockIdx++) {
        const LLKA_CifDataBlock *block = &data->blocks[blockIdx];
        const LLKA_CifDataCategory *cat = block->firstCategory;

        fprintf(stdout, "Data frame %s\n", block->name);

        while (cat) {
            const LLKA_CifDataItem *item = cat->firstItem;

            fprintf(stdout, "_%s\n", cat->name);

            while (item) {
                size_t valIdx;

                fputc('\t', stdout);
                fputs(item->keyword, stdout);
                for (valIdx = 0; valIdx < item->nValues; valIdx++) {
                    const LLKA_CifDataValue *v = &item->values[valIdx];

                    putchar(' ');
                    if (v->state == LLKA_MINICIF_VALUE_SET)
                        fputs(v->text, stdout);
                    else if (v->state == LLKA_MINICIF_VALUE_NONE)
                        fputc('.', stdout);
                    else if (v->state == LLKA_MINICIF_VALUE_UNKW)
                        fputc('?', stdout);
                }
                fputc('\n', stdout);

                item = LLKA_cifDataCategory_nextItem(item);
            }

            fputc('\n', stdout);
            cat = LLKA_cifDataBlock_nextCategory(cat);
        }

        fputc('\n', stdout);
    }
}

int main(int argc, char *argv[])
{
    LLKA_CifData *data;
    LLKA_RetCode tRet;
    char *cifError;

    if (argc < 2) {
        fprintf(stderr, "No input file\n");

        return EXIT_FAILURE;
    }

    tRet = LLKA_cifFileToData(argv[1], &data, &cifError);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to read Cif file: %s\n", LLKA_errorToString(tRet));
        if (cifError) {
            fprintf(stderr, "%s\n", cifError);
            LLKA_destroyString(cifError);
        }

        return EXIT_FAILURE;
    }

    mirror(data);

    LLKA_destroyCifData(data);

    return EXIT_SUCCESS;
}
