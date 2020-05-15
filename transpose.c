#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MEM_OVERHEAD 1.01    /* speed hack -- overallocate to avoid reallocs */

struct text_data
{
    char     **lines;
    void     **line_idxs;
    uint8_t   *line_bits;
    uint32_t  *num_cols_per_row;
    uint32_t  *num_rows_per_col;
    uint32_t   num_rows;
    uint32_t   num_cols;
};


/* realloc input string, return new max array length (including terminal null)
 * handles \r\n \n \r, including mixes of EOL characters within same file
 * strips EOL from end of string
 */
char * fgets_strip_realloc_len(char **return_string,
                               uint32_t *return_max_length, FILE *infile,
                               uint32_t *return_length)
{
    int c;
    char     *string        = *return_string;
    uint32_t  length        = 0;
    uint32_t  total_length;
    uint32_t  max_length    = *return_max_length;
    int       old_c         = '\0';
    uint32_t  anything_flag = 0;

    while((c = fgetc(infile)) != EOF)
    {
        anything_flag = 1;
    
        /* EOL: \n or \r\n */
        if (c == '\n')
        {
            /* MSDOS, get rid of the previously stored \r */
            if (old_c == '\r')
            {
                string[length-1] = '\0';
            }

            old_c = c;
            
            break;
        }
        /* EOL: \r */
        /* may be a Mac text line, back up a character */
        else if (old_c == '\r')
        {
            /* fseek(infile, -1 * sizeof(char), SEEK_CUR); */
            ungetc(c, infile);

            break;
        }
        
        old_c = c;
    
        total_length = length + 2;    /* increment, plus terminal null */

        if (total_length > max_length)
        {
            max_length = MEM_OVERHEAD * total_length;
            string     = (char *) realloc(string, max_length * sizeof(char));
        }

        string[length++] = c;
    }
    
    /* check for dangling \r from reading in Mac lines */
    if (length && string[length-1] == '\r')
    {
        length--;
    }
    
    if (length == 0)
    {
        if (!anything_flag)
        {
            *return_length = 0;

            return NULL;
        }
            
        if (max_length < 1)
        {
            max_length = 1;
            string     = (char *) realloc(string, sizeof(char));
        }
    }

    string[length]     = '\0';
    
    *return_string     = string;
    *return_max_length = max_length;
    *return_length     = length;

    return string;
}


/* replace tabs with nulls, fill array of field pointers,
 * return number of fields
 * WARNING -- clobbers tabs in original input string
 */
uint32_t split_tabs(char *string, char ***fields, uint32_t *return_max_field)
{
    char *cptr, *sptr;
    uint32_t count     = 0;
    uint32_t max_field = *return_max_field;
    
    sptr = string;
    for (cptr = string; *cptr; cptr++)
    {
        if (*cptr == '\t')
        {
            count++;

            if (count > max_field)
            {
                max_field = MEM_OVERHEAD * count;
                *fields   = realloc(*fields, max_field * sizeof(char *));
            }

            (*fields)[count-1] = sptr;
            sptr = cptr + 1;
            *cptr = '\0';
        }
    }
    
    /* final field */
    count++;
    if (count > max_field)
    {
        max_field = MEM_OVERHEAD * count;
        *fields = realloc(*fields, max_field * sizeof(char *));
    }
    (*fields)[count-1] = sptr;
    
    *return_max_field = max_field;
    
    return count;
}


int read_text_data(char *filename, struct text_data *text_data)
{
    FILE *infile;
    uint32_t max_string_len = 0;
    char *string            = NULL;
    char **fields           = NULL;
    uint32_t num_fields     = 0;
    uint32_t max_num_fields = 0;
    char *buffer            = NULL;
    uint32_t i, k;
    uint32_t row;
    
    char     **lines                = NULL;
    void     **line_idxs            = NULL;
    uint32_t  *ptr_u32;
    uint16_t  *ptr_u16, u16;
    uint8_t   *ptr_u8;
    uint8_t   *line_bits            = NULL;
    uint32_t  *num_cols_per_row     = NULL;
    uint32_t  *num_rows_per_col     = NULL;
    uint32_t   max_allocated_rows   = 0;
    uint32_t   max_allocated_cols   = 0;
    uint32_t   old_allocated_cols   = 0;
    uint32_t   max_row = 0, max_col = 0;
    uint32_t   line_length, max_field_index;
    
    /* allocate my own i/o buffer, since we can't trust the system and/or
     * compiler to allocate a decently large one...
     */
    
    buffer = (char *) malloc(1048576 * sizeof(char));

    /* standard input */
    if (filename == NULL || strcmp(filename, "-") == 0)
    {
        infile = stdin;
    }
    else
    {
        infile = fopen(filename, "rb");
        if (!infile)
        {
            fprintf(stderr, "ERROR -- can't open input file %s\n", filename);

            if (buffer) free(buffer);
            if (string) free(string);
            if (fields) free(fields);
            
            return 1;
        }
        setvbuf(infile, buffer, _IOFBF, 1048576);
    }

    row = 0;
    while(fgets_strip_realloc_len(&string, &max_string_len, infile, &line_length))
    {
        num_fields = split_tabs(string, &fields, &max_num_fields);

        /* allocate more col row counters, 1% overhead for efficiency */
        if (num_fields > max_allocated_cols)
        {
            old_allocated_cols = max_allocated_cols;
            max_allocated_cols = MEM_OVERHEAD * num_fields;
            num_rows_per_col = (uint32_t *) realloc(num_rows_per_col,
                                        max_allocated_cols * sizeof(uint32_t));
            
            /* zero out the newly allocated counters */
            memset(&num_rows_per_col[old_allocated_cols], 0,
                   (max_allocated_cols - old_allocated_cols) *
                   sizeof(uint32_t));
        }

        /* allocate more row col counters, 1% overhead for efficiency */
        if (row + 1 > max_allocated_rows)
        {
            max_allocated_rows = MEM_OVERHEAD * (row + 1);

            lines              = (char **)    realloc(lines,
                                                      max_allocated_rows *
                                                      sizeof(char *));
            line_idxs          = (void **)    realloc(line_idxs,
                                                      max_allocated_rows *
                                                      sizeof(void *));
            num_cols_per_row   = (uint32_t *) realloc(num_cols_per_row,
                                                      max_allocated_rows *
                                                      sizeof(uint32_t));
            line_bits          = (uint8_t *)  realloc(line_bits,
                                                      max_allocated_rows *
                                                      sizeof(uint8_t));
        }

        /* check for no real fields at all */
        line_idxs[row] = NULL;
        if (num_fields == 1 && fields[0][0] == '\0')
            num_fields = 0;

        num_cols_per_row[row] = num_fields;

        
        /* I have to be clever here, otherwise memory allocation overhead
         *  explodes the memory usage.
         * Store the substrings as one big null-delimited line,
         *  and store separate indices to the substrings within it.
         * This way, we only do two malloc() per line, rather than
         *  a whole lot.  Memory savings is *SUBSTANTIAL* on big files.
         *
         * #rows and #cols are limited to <= 2^32 each, due to using
         *  32-bit unsigned integer book keeping.
         *
         * Maximum line length per row is also limited to <= 2^32, in order
         *  to use 32-bit field indices, rather than 64-bit field pointers.
         * Line input functions were already using 32-bit integers, so we
         *  aren't actually losing anything vs. using pointer arrays instead.
         *
         * I'm guessing from my examples so far that worst-case memory usage
         *  will wind up being no more than ~1.5x the original file size for
         *  typical large files containing mostly numerical data.  Very small
         *  files (tens of megabytes) will have worse ratios, due to the
         *  book keeping arrays being of similar size to the data itself.
         */
        if (num_fields)
        {
            lines[row] =
                (char *) malloc((line_length + 1) * sizeof(char));
            memcpy(lines[row], string, (line_length + 1) * sizeof(char));


            max_field_index = fields[num_fields-1] - string;
        
            /* standard 32-bit indices */
            if (max_field_index > 0xFFFFFF)
            {
                line_idxs[row] =
                    (uint32_t *) malloc(num_fields * sizeof(uint32_t));

                /* store fields in row */
                ptr_u32 = (uint32_t *) line_idxs[row];
                for (i = 0; i < num_fields; i++)
                    ptr_u32[i] = fields[i] - string;
                
                line_bits[row] = 32;
            }
            /* HACK -- cram into 24-bit indices to save memory */
            else if (max_field_index > 0xFFFF)
            {
                line_idxs[row] =
                    (uint8_t *) malloc(num_fields * 3*sizeof(uint8_t));

                /* store fields in row */
                ptr_u8 = (uint8_t *) line_idxs[row];
                for (i = 0; i < num_fields; ptr_u8 += 3, i++)
                {
                    k         = fields[i] - string;
                    u16       = (uint16_t) k;
                    memcpy(ptr_u8, &u16, sizeof(uint16_t));
                    ptr_u8[2] = k >> 16;
                }
                
                line_bits[row] = 24;
            }
            /* HACK -- cram into 16-bit indices to save memory */
            else if (max_field_index > 0xFF)
            {
                line_idxs[row] =
                    (uint16_t *) malloc(num_fields * sizeof(uint16_t));

                /* store fields in row */
                ptr_u16 = (uint16_t *) line_idxs[row];
                for (i = 0; i < num_fields; i++)
                    ptr_u16[i] = fields[i] - string;
                
                line_bits[row] = 16;
            }
            /* HACK -- cram into 8-bit indices to save memory */
            else
            {
                line_idxs[row] =
                    (uint8_t *) malloc(num_fields * sizeof(uint8_t));

                /* store fields in row */
                ptr_u8 = (uint8_t *) line_idxs[row];
                for (i = 0; i < num_fields; i++)
                    ptr_u8[i] = fields[i] - string;
                
                line_bits[row] = 8;
            }
        }
        else
        {
            line_idxs[row] = NULL;
            lines[row]     = NULL;
        }

        row++;

        /* store number of rows for each col */
        num_rows_per_col[0] = row;
        for (i = 1; i < num_fields; i++)
            num_rows_per_col[i] = row;
            
        if (row > max_row)
            max_row = row;
        if (num_fields > max_col)
            max_col = num_fields;
    }
    
    fclose(infile);
    
    text_data->lines            = lines;
    text_data->line_idxs        = line_idxs;
    text_data->line_bits        = line_bits;
    text_data->num_cols_per_row = num_cols_per_row;
    text_data->num_rows_per_col = num_rows_per_col;
    text_data->num_rows         = max_row;
    text_data->num_cols         = max_col;

    if (buffer) free(buffer);
    if (string) free(string);
    if (fields) free(fields);
    
    return 0;
}


void transpose_text_data(struct text_data *text_data, FILE *outfile)
{
    char     **lines            = text_data->lines;
    void     **line_idxs        = text_data->line_idxs;
    char      *string;
    uint32_t  *ptr_u32;
    uint16_t  *ptr_u16, u16;
    uint8_t   *ptr_u8;
    uint8_t   *line_bits        = text_data->line_bits;
    uint32_t  *num_cols_per_row = text_data->num_cols_per_row;
    uint32_t  *num_rows_per_col = text_data->num_rows_per_col;
    uint32_t   num_rows         = text_data->num_rows;
    uint32_t   num_cols         = text_data->num_cols;
    uint32_t   row, col;

    /* file is empty, abort */
    if (num_rows == 0)
        return;
        
    /* print first new row, since we don't need to look up field indices */
    for (row = 0; row < num_rows_per_col[0]; row++)
    {
        if (row)
            fputc('\t', outfile);

        if (num_cols_per_row[row])
            fputs(lines[row], outfile);
    }
    fputc('\n', outfile);

    
    /* print the rest of the new rows */
    for (col = 1; col < num_cols; col++)
    {
        for (row = 0; row < num_rows_per_col[col]; row++)
        {
            if (row)
                fputc('\t', outfile);
            
            /* print data present for col in this row */
            if (col < num_cols_per_row[row])
            {
                /* NOTE -- converting into a switch/case was slower (!!) */
            
                /* HACK -- cram into 24-bit indices to save memory */
                /* most common case for files large enough to matter */
                if (line_bits[row] == 24)
                {
                    ptr_u8  = (uint8_t *) line_idxs[row] + 3*col;
                    memcpy(&u16, ptr_u8, sizeof(uint16_t));
                    string  = lines[row] +
                              ((uint32_t) u16 | (ptr_u8[2] << 16));
                }
                /* HACK -- cram into 16-bit indices to save memory */
                /* 2nd most common case */
                else if (line_bits[row] == 16)
                {
                    ptr_u16 = (uint16_t *) line_idxs[row];
                    string  = lines[row] + ptr_u16[col];
                }
                /* HACK -- cram into 8-bit indices to save memory */
                /* distant 3rd most common case, so fast it won't matter */
                else if (line_bits[row] == 8)
                {
                    ptr_u8  = (uint8_t *) line_idxs[row];
                    string  = lines[row] + ptr_u8[col];
                }
                /* standard 32-bit indices */
                /* we don't expect to encounter lines this long */
                else
                {
                    ptr_u32 = (uint32_t *) line_idxs[row];
                    string  = lines[row] + ptr_u32[col];
                }

                fputs(string, outfile);
            }
        }
        fputc('\n', outfile);
    }
}


void free_text_data(struct text_data *text_data)
{
    char     **lines            = text_data->lines;
    void     **line_idxs        = text_data->line_idxs;
    uint8_t  *line_bits         = text_data->line_bits;
    uint32_t  *num_cols_per_row = text_data->num_cols_per_row;
    uint32_t  *num_rows_per_col = text_data->num_rows_per_col;
    uint32_t   num_rows         = text_data->num_rows;
    uint32_t   row;
    
    if (line_idxs)
    {
        for (row = 0; row < num_rows; row++)
        {
            if (lines[row])
                free(lines[row]);

            if (line_idxs[row])
                free(line_idxs[row]);
        }

        if (lines)
            free(lines);

        if (line_idxs)
            free(line_idxs);
    }

    if (line_bits)        free(line_bits);
    if (num_cols_per_row) free(num_cols_per_row);
    if (num_rows_per_col) free(num_rows_per_col);
}


int transpose_entry_point(char *infile_name, char *outfile_name)
{
    FILE   *outfile    = NULL;
    char   *buffer_out = NULL;
    struct  text_data text_data;
    
    if (read_text_data(infile_name, &text_data))
    {
        return(1);
    }

    /* open file for writing */    
    if (outfile_name)
    {
        /* open as text, so it will translate EOL automatically */
        outfile = fopen(outfile_name, "wt");

        if (!outfile)
        {
            fprintf(stderr, "ERROR -- can't open output file %s\n",
                    outfile_name);

            free_text_data(&text_data);

            return 1;
        }
    }
    /* write to stdout */
    else
    {
        outfile = stdout;
    }

    buffer_out = (char *) malloc(1048576 * sizeof(char));
    setvbuf(outfile, buffer_out, _IOFBF, 1048576);

    transpose_text_data(&text_data, outfile);
    fclose(outfile);
    
    free(buffer_out);
    free_text_data(&text_data);
    
    return 0;
}

int main (int argc, char *argv[])
{
    char *infile_name  = NULL;
    char *outfile_name = NULL;

    if (argc > 1)
    {
        infile_name  = argv[1];
    }
    if (argc > 2)
    {
        outfile_name = argv[2];
    }
    
    if (transpose_entry_point(infile_name, outfile_name))
    {
        /* ERROR -- exit with non-zero value */
        return 1;
    }
    
    return 0;
}
