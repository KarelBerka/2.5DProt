#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <cjson/cJSON.h>

#define M_PI 3.14159265358979323846

typedef struct {
    char *name;
    void *pointer;
} StrList;

typedef struct {
   char *name;
   float x;
   float y;
   void *pointer;
} XYLayout;

typedef struct {
   char *name;
   float angle;
   void *pointer;
} AngleLayout;

typedef struct {
   char *name;
   float x;
   float y;
   float z;
   void *pointer;
} ThreeDLayout;

typedef struct {
   char *name;
   float diff_value;
   float twod_diff;
   float size;
   void *pointer;
} Diff;

typedef struct {
  char *name;
  void *row;
  Diff *pointer;
} ChDiff;

typedef struct {
   char *name;
   float row_size;
   void *row;
   ChDiff *pointer;
} TabDiff;

// aux functions


char *getSubString(int start, int end, char string[])
{
    char *substring = (char *)malloc(sizeof(char ) * (end -start +1));
    int len = end - start;

    for (int i=0;i < len;i++) {
        substring[i] = string[start+i];
    }
    substring[len] = '\0';
    return substring;
}

// StrList functions

char *strlist_find_nth_value(int n, StrList *values)
{
    StrList *rec = values;
    int counter = 0;

    while (counter < n) {
         rec = rec->pointer;
         counter = counter +1;
    }
    if (rec == NULL) {
        printf("ERROR: strlist_find_nth_value %d\n", n);
    }
    return rec->name;
}


int strlist_length(StrList *values)
{
    StrList *rec = values;
    int counter = 0;

    while (rec) {
         rec = rec->pointer;
         counter = counter +1;
    }
    return counter;
}


// print structures

char print_list(char *name, StrList *values)
{
    StrList *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s, ", rec->name);
        rec = rec->pointer;
    }
    printf("\n");
}


char print_xylayout(char *name, XYLayout* values)
{
    XYLayout *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s<x:%.2f,y:%.2f>, ", rec->name, rec->x, rec->y);
        rec = rec->pointer;
    }
    printf("\n");
}


char print_anglelayout(char *name, AngleLayout* values)
{
    AngleLayout *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s<angle:%.2f>, ", rec->name, rec->angle);
        rec = rec->pointer;
    }
    printf("\n");
}


char print_3dlayout(char *name, ThreeDLayout* values)
{
    ThreeDLayout *rec = values;
    printf("\n-- %s --: ", name);
    while (rec != NULL) {
        printf("%s<x:%4.1f,y:%4.1f,z:%5.1f>, ", rec->name, rec->x, rec->y, rec->z);
        rec = rec->pointer;
    }
    printf("\n");
}


char print_tabdiff(char *name, TabDiff* values)
{
    // write table - each row contains data from one row instance
    // all its diff values 
    int c = 0;
    
    Diff *diff, *header;
    ChDiff *chdiff;
    TabDiff *row;
    printf("-- %s --: \n", name);
    row = values;
    while (row != NULL) {
        chdiff = row->pointer;
        printf("%4s ", row->name);
        while (chdiff != NULL) {
            diff = chdiff->pointer;
            while (diff != NULL) {
                printf("{%2.0f}x{%2.0f}(%s) ", diff->diff_value, diff->twod_diff, diff->name);
                diff = diff ->pointer;
            }
            chdiff = chdiff->row;
        }
        printf("konec \n");
        row = row->row;
    }
}

char print_tabdiff_2(char *name, TabDiff* values)
{
    // write table - each row contains data from one row instance
    // all its diff values 
    int c = 0;

    Diff *diff, *header;
    ChDiff *chdiff;
    TabDiff *row;
    printf("-- %s --: \n", name);
    row = values;
    while (row != NULL) {
        chdiff = row->pointer;
        printf("===========================\n");
        printf("  %4s ", row->name);
        while (chdiff != NULL) {
            diff = chdiff->pointer;
            printf("   onelevel: %s -", chdiff->name);
            while (diff != NULL) {
                printf("{%2.0f}x{%2.0f}:{%2.0f} (%s), ", diff->diff_value, diff->twod_diff, diff->size, diff->name);
                diff = diff ->pointer;
            }
            chdiff = chdiff->row;
            printf("\n");
        }
        row = row->row;
    }
}




// --------------------------------------------------------------------------------------------------------------------------------------
// read input structures
// --------------------------------------------------------------------------------------------------------------------------------------
StrList *read_str_list_variable(cJSON *cjson_buffer_var)
{
    cJSON *child;
    StrList *list = NULL;
    StrList *rec = NULL, *new_rec;

    child = cjson_buffer_var->child;
    while (child != NULL) {
        new_rec = (StrList *)malloc((sizeof(StrList)+1));   // create new record rec->name is the name, pointer to null
        new_rec->name = child->valuestring;
        new_rec->pointer = NULL;
        if (list == NULL) {                             // the first record, mark it as start
            list = new_rec;
        } else {                                        // not the first, connect it to the line
            rec->pointer = new_rec;
        }
        rec = new_rec;
        child = child->next;
    }
    return list;
}


XYLayout *read_xylayout_variable(cJSON *cjson_buffer_var)
{
    XYLayout *list = NULL;
    XYLayout *rec = NULL, *new_rec;

    cJSON *child, *coords;
    child = cjson_buffer_var->child;

    while (child != NULL) {
        coords = child->child;
        new_rec = (XYLayout *)malloc(sizeof(XYLayout)); // create new record rec->name is the name, pointer to null
        new_rec->name = child->string;
        new_rec->x = coords->valuedouble;
        coords = coords-> next;
        new_rec->y = coords->valuedouble;
        new_rec->pointer = NULL;
        if (list == NULL) {                             // the first record, mark it as start
            list = new_rec;
        } else {                                        // not the first, connect it to the line
            rec->pointer = new_rec;
        }
        rec = new_rec;
        child = child->next;
    }
    return list;
}


AngleLayout *read_anglelayout_variable(cJSON *cjson_buffer_var)
{
    AngleLayout *list = NULL;
    AngleLayout *rec = NULL, *new_rec;

    cJSON *child;
    child = cjson_buffer_var->child;
    int val;

    while (child != NULL) {
        new_rec = (AngleLayout *)malloc(sizeof(AngleLayout)); // create new record rec->name is the name, pointer to null
        new_rec->name = child->string;
        new_rec->angle = (float)(child->valueint);
        new_rec->pointer = NULL;
        if (list == NULL) {                             // the first record, mark it as start
            list = new_rec;
        } else {                                        // not the first, connect it to the line
            rec->pointer = new_rec;
        }
        rec = new_rec;
        child = child->next;
    }
    return list;
}


ThreeDLayout *read_3dlayout_variable(char *LineBuffer, char *name)
{
    char *pos;
    int sln = strlen(name);

    ThreeDLayout *list = NULL;
    ThreeDLayout *rec = NULL, *new_rec;

    char *sse_name;
    int sse_start = 0, sse_end = 0, val_start = 0, val_middle = 0, val_end;

    if (strlen(LineBuffer) <= sln) {
        return list;
    }
    if (strncmp(LineBuffer, name, sln) != 0) {
        return list;
    }

    pos = strchr(LineBuffer, ' ');
    for (int i=0; i<strlen(pos); i++) {
        if (pos[i] == '\'') {                             // sse name is read (start/end)
           if (sse_start == 0) {
               sse_start = i+1;
           } else {
               sse_end = i;
           }
        }

        if ((pos[i] == '(') || (pos[i] == '[')) {                              // value is read - start
            val_start = i+1;
        }
        if ((pos[i] == ',') && (val_start != 0)) {        // value is read - middle part
            val_middle = i;
        }
        if ((pos[i] == ',') && (val_end != 0)) {        // value is read - middle part
            val_end = i;
        }
        if ((pos[i] == ')') || (pos[i] == ']')){                              // value is read - end
            new_rec = (ThreeDLayout *)malloc(sizeof(XYLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = getSubString(sse_start, sse_end, pos);
            new_rec->x = atof(getSubString(val_start, val_middle, pos));
            new_rec->y = atof(getSubString(val_middle+1, val_end, pos));
            new_rec->z = atof(getSubString(val_end+1, i, pos));
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
            } else {                                        // not the first, connect it to the line
                rec->pointer = new_rec;
            }
            rec = new_rec;
            sse_end = 0;
            sse_start = 0;
            val_start = 0;
            val_middle = 0;
            val_end = 0;
        }
    }
    return list;
}


char *sses_channel_prefix(char *name)
{
    int counter = 0;
//    char ret[strlen(name)+1];
    while (counter<strlen(name)) {
        if (isalpha(name[counter])) {
            if ((counter+1<strlen(name)) && (isdigit(name[strlen(name)-1]))) {            
//                strcpy(ret, name);
//                printf("ret je %s konci na: %c, z toho je prefix ->", ret, ret[strlen(ret)-1]);
                return getSubString(0, counter, name);
//                printf("prefix je %s\n", ret);
            } else {
                return "";
            }
        }
        counter += 1;
    }
    return "";
}


TabDiff *read_tabdiff_variable(cJSON *cjson_buffer_var)
{
    TabDiff *row_start = NULL, *row = NULL , *row_new = NULL;
    Diff *diff = NULL, *diff_new = NULL, *cmp_pointer = NULL;
    ChDiff *chdiff, *chdiff_new;

    cJSON *cjson_row, *child;
    cjson_row = cjson_buffer_var->child;
    int res;
    int leave;
    char *prefix;
    char *ret;

    while (cjson_row != NULL) {
        child = cjson_row->child;
        row_new = (TabDiff *)malloc(sizeof(TabDiff));        // create new row instance (name, null, null)
        row_new->name = cjson_row->string;
        row_new->row = NULL;
        row_new->pointer = NULL;
        if (row_start == NULL) {                             // there was no previous - set is at start
            row_start = row_new;
        } else {
            row->row = row_new;                              // there was previous - connect the new one to old->row
        }
        row = row_new;
        chdiff = NULL;
        while (child != NULL) {
            prefix = sses_channel_prefix(child->string);
            // there is no need to go through whole chdiff list channels are sorted
            if ((chdiff == NULL) || (strcmp(chdiff->name, prefix) != 0)) {
                // we have new channel line
                chdiff_new = (ChDiff *)malloc(sizeof(ChDiff));
                chdiff_new->name = prefix;
                chdiff_new->pointer = NULL;
                chdiff_new->row = NULL;
                if (row->pointer == NULL) {
                    row->pointer = chdiff_new;
                } else {
                    chdiff->row = chdiff_new;
                }
                chdiff = chdiff_new;
            }

            diff_new = (Diff *)malloc(sizeof(Diff));  // create new diff instance (name,value, null, null)
            diff_new->name = child->string;
            diff_new->diff_value = child->valueint;
            diff_new->pointer = NULL;
            if (chdiff->pointer == NULL) {
                // start of the new row for new sse
                chdiff->pointer = diff_new;
            } else {
                // we have to find proper position
                diff->pointer = diff_new;
            }
            diff = diff_new;
            // diff record is on the right place in the tree
            child = child->next;
        }
        cjson_row = cjson_row->next;
    }
    return row_start;
}


// --------------------------------------------------------------------------------------------------------------------------------------
// all data read
// --------------------------------------------------------------------------------------------------------------------------------------
/*
float find_sse_coords_firstonly(TabDiff *lsdistances, char *chan_name, char *sse_name)
{
    TabDiff *row = lsdistances;
    Diff *cell;
    while (row) {
        if (strcmp(row->name, chan_name) == 0) {
            cell = row->pointer;
            while (cell) {
                if  (strcmp(cell->name, sse_name) == 0) {
                   return cell->dif_value;
                }
                cell = cell->pointer;
            }
        }
        row = row->row;
    }
    return -1;
}
*/
TabDiff *find_sse_coords_row(TabDiff *lsdistances, char *chan_name)
{
    TabDiff *row = lsdistances;
    while (row) {
        if (strcmp(row->name, chan_name) == 0) {
            return row;
        }
        row = row->row;
    }
    return NULL;
}

float find_sse_coords_from_row(TabDiff *row, char *sse_name)
{
    ChDiff *ch_row = row->pointer;
    Diff *cell;
    char *prefix;

    prefix = sses_channel_prefix(sse_name);

    while (ch_row) {

        if (strcmp(ch_row->name, prefix) == 0) {
            break;
        }
        ch_row = ch_row->row;
    }

    if (strcmp(ch_row->name, prefix) != 0) {
        return -1;
    }
    cell = ch_row->pointer;
    while (cell) {
        if  (strcmp(cell->name, sse_name) == 0) {
            return cell->diff_value;
        }
        cell = cell->pointer;
    }
    return -1;
}


int is_string_in_strlist(char *name, StrList *list)
{
    StrList *rec = list;
    while(rec != NULL) {
        if (strcmp(rec->name, name) == 0) {
            return(0);
        }
        rec = rec->pointer;
    }
    return(1);
}


StrList *char_to_list(char *for_helic)
{
    StrList *rec_new, *rec = NULL, *list = NULL;
    int counter = 0;
    int start = 0;

    while (counter < strlen(for_helic)) {
        if (for_helic[counter] == ':') {    // we have sse name
            rec_new = (StrList *)malloc(sizeof(StrList));
            rec_new->name = getSubString(start, counter, for_helic);
            rec_new->pointer = NULL;
            if (list == NULL) {
                list = rec_new;
            } else {
                rec->pointer = rec_new;                              // row only continue
            }
            rec = rec_new;
            start = counter + 1;
        }
        counter = counter + 1;
    }
    rec_new = (StrList *)malloc(sizeof(StrList));
    rec_new->name = getSubString(start, counter, for_helic);
    rec_new->pointer = NULL;
    if (list == NULL) {
        list = rec_new;
    } else {
        rec->pointer = rec_new;                              // row only continue
    }
    rec = rec_new;

    return list;
}


int change_xylayout_with_diffs(XYLayout *full_layout, char *name, float dx, float dy)
{
     XYLayout *res = full_layout;
     
     while (res) {
         if (strcmp(name, res->name) == 0) {
//             printf("pred %s %f %f po %f %f ", res->name, res->x, res->y, dx, dy);
             res->x = res->x +dx;
             res->y = res->y +dy;
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}


int print_xylayout_value(XYLayout *full_layout, char *name)
{
     XYLayout *res = full_layout;
     
     while (res) {
         if (strcmp(name, res->name) == 0) {
             printf(" (%f, %f)", res->name, res->x, res->y);
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}




float count_channel_error(StrList *for_channel_list, XYLayout *channels_layout, StrList *all_s, StrList *prime_s, AngleLayout *sses_sizes, XYLayout *sses_layout, TabDiff *lsdistances)
{

    StrList *rec=for_channel_list;
    XYLayout *r2=channels_layout;
    XYLayout *sse=sses_layout;
    float l_dist, td_dist, min_dist;
    float d_diff, lt_diff, diff = 0;
    float error = 0;
    int ccc = 0;
    XYLayout *ch_layout;
    TabDiff *row;

    while (rec) {
        ch_layout = channels_layout;
        row = find_sse_coords_row(lsdistances, rec->name);
        if (row == NULL) {
             printf("name %s was not found %p\n", rec->name, row);
             continue;
        }
        while (ch_layout) {
            if (strcmp(ch_layout->name, rec->name) == 0) {
                break;
            }
            ch_layout = ch_layout->pointer;
        }
//        printf("radek ./ ", rec->name);
        
        while (sse) {
            l_dist = sqrt(((ch_layout->x - sse->x)*(ch_layout->x - sse->x))+ ((ch_layout->y - sse->y)*(ch_layout->y - sse->y)));
//            td_dist =find_sse_coords_firstonly(lsdistances, rec->name, sse->name);
//            printf("pred fscfr\n");
            td_dist = find_sse_coords_from_row(row, sse->name);
//            printf("po fscrf %f %s %s\n", td_dist, sse->name, rec->name);
            min_dist = fmin(l_dist, td_dist);
            lt_diff = (fabs(l_dist-td_dist) + 3 * fmax(0, l_dist-td_dist))/4;
            d_diff = 0;

            if (min_dist < 25) {
                d_diff = lt_diff * (9 - min_dist/4) * (9 - min_dist/4);
            }

            if ((min_dist >= 25) && (min_dist < 35)) {
                d_diff = lt_diff / 2;
            }

            if ((strlist_length(prime_s) >= 2) && (!(is_string_in_strlist(sse->name, prime_s)))) {
                d_diff = d_diff / 2;
            }
//            printf("diff(%d x %d->%d)", (int) l_dist, (int) td_dist,(int)d_diff);
            ccc += 1;
            diff += d_diff;
            sse = sse->pointer;
        }
//        printf("sses error %f (%d) - ", diff, ccc);
        while (r2) {
            l_dist = sqrt(((ch_layout->x - r2->x)*(ch_layout->x - r2->x))+ ((ch_layout->y - r2->y)*(ch_layout->y - r2->y)));
            td_dist = find_sse_coords_from_row(row, r2->name);
//            td_dist = 3;
            min_dist = fmin(l_dist, td_dist);
//            lt_diff = (fabs(l_dist-td_dist) + 3 * fmax(0, l_dist-td_dist))/4;
///            lt_diff = fmax(0, l_dist-td_dist);
            if (l_dist < td_dist) {
                lt_diff = (td_dist -l_dist)/4;
            } else {
               lt_diff = l_dist-td_dist;
            }
            if (min_dist < 10) {
                diff += lt_diff * 3;
                d_diff = lt_diff * 3;
            }
            if ((min_dist >= 10) && (min_dist < 25)) {
                diff += lt_diff;
                d_diff = lt_diff;
            }
            r2 = r2->pointer;
            ccc += 1;
//            p1 = (int)lt_diff;
//            p2 = (int)d_diff;
//            printf("diff(%d x %d->%d)", (int) l_dist, (int) td_dist,(int)d_diff);
        }
        rec = rec->pointer;
//        printf("full error %f (%d)\n", diff, ccc);

    }
    return diff;
}


int iteration(StrList *ligandsnames, XYLayout *channels_layout, StrList *all_s, StrList *prime_s, AngleLayout *sses_sizes, XYLayout *sses_layout, TabDiff *lsdistances, int STEP_SIZE, int CHECK_EVERY)
{
    int max_generations = 500000;
//    int max_generations = 10000;
    float MINIMAL_IMPROVEMENT = 0.03;   // if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout

    long int generation = 1;
    float error_checkpoint = -1;
    int coef = 1;

    float error_old, error_new, wrror_checkpoint = 0;

    int len_ligandsnames = strlist_length(ligandsnames);
    int r;
    char *for_channel;
    StrList *for_channel_list;
    float dx, dy;
    float errors_summary, errors_summary_old;

    float cxmin, cxmax, cymin, cymax;
    float eavg, esum = 0;
    AngleLayout *aux, *err_start = NULL, *err_end = NULL;
    XYLayout *lname;
    XYLayout *ch_layout_fix;
    int counter;

    // we want to have the stable results for all iterations
//    srand(time(NULL));                      // initialize random function


    errors_summary_old = count_channel_error(ligandsnames, channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);
    while (generation <= max_generations+1) {
        // choose the moved ligand
        r = rand() % len_ligandsnames;
        for_channel = strlist_find_nth_value(r, ligandsnames);
        for_channel_list = char_to_list(for_channel);

        error_old = count_channel_error(for_channel_list, channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);
        // count new error
        dx = ((float)(rand() % ((int)(STEP_SIZE/coef*10))))/10 - STEP_SIZE/(2*coef);
        dy = ((float)(rand() % ((int)(STEP_SIZE/coef*10))))/10 - STEP_SIZE/(2*coef);
        change_xylayout_with_diffs(channels_layout, for_channel, dx, dy);

        error_new = count_channel_error(for_channel_list, channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);

        if (error_new >= error_old) {
            // return coordinates to new one
            change_xylayout_with_diffs(channels_layout, for_channel, -dx, -dy);

        } else {
            printf("%ld: LIGANDS change %s o dx:(%f) dy:(%f) max(%d) err before(%f) after(%f)\n", generation, for_channel, dx, dy, (STEP_SIZE/coef)/2, error_old, error_new);
//            print_xylayout_value(channels_layout, for_channel);
            // set the error to new one
            error_old = error_new;
        }
        if ((generation % CHECK_EVERY) == 0) {
            errors_summary = count_channel_error(ligandsnames, channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);
//            if (errors_summary != 0) {
//                printf("%d stara chyba: %f nova: %f procenta: %f - che(%d)\n", generation, errors_summary_old, errors_summary, ((errors_summary_old - errors_summary)/errors_summary*100), CHECK_EVERY);
//            } else {
//                printf("error 0\n");
//            }  

            if (generation == CHECK_EVERY) {
//                 If we are in the first round, there might be "lost" ligands which are in absolutely different position
//                 if there is a lot of ligands, we have to put them to the best positions

//                 count average error of all channels
                   esum = 0;
                   cxmin = 500;
                   cxmax = 0;
                   cymin = 500;
                   cymax = 0;
                   lname = channels_layout;
                   while (lname) {
                       aux = (AngleLayout *)malloc(sizeof(AngleLayout)+1);
                       aux->name = lname->name;
                       aux->angle = count_channel_error(char_to_list(lname->name), channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);
                       aux->pointer = NULL;
                       if (err_start == NULL) {
                           err_start = aux;
                           err_end = aux;
                       } else {
                           err_end->pointer = aux;
                       }
                       esum = esum + aux->angle;
                       cxmin = fmin(cxmin, lname->x);
                       cxmax = fmax(cxmax, lname->x);
                       cymin = fmin(cymin, lname->y);
                       cymax = fmax(cymax, lname->y);
                       lname = lname->pointer;
                   }
                   eavg = esum / len_ligandsnames;
                   aux = err_start;
                   while (aux) {
                       if (aux->angle > eavg * 1.8) {
                           // if the error is bigger then twice the standart
                           // - we the ligand could have wrong position
                           ch_layout_fix = channels_layout;
                           while (ch_layout_fix) {
                               if (strcmp(ch_layout_fix->name,aux->name) == 0) {
                                   // tru to find better layout for this channel position
                                   dx = ch_layout_fix->x;
                                   dy = ch_layout_fix->y;
                                   for (r = 0; r < 100; r ++) {
                                       //go through 10 to 10 and try them
                                       ch_layout_fix->x = cxmin + (cxmax - cxmin) * (r % 10) / 9;
                                       ch_layout_fix->y = cymin + (cymax - cymin) * ((int) r / 10) / 9;
                                       for_channel_list = char_to_list(aux->name);
                                       error_new = count_channel_error(for_channel_list, channels_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances);
                                       if (error_new < aux->angle) {
                                           aux->angle = error_new;
                                       } else {
                                           ch_layout_fix->x = dx;
                                           ch_layout_fix->y = dy;
                                       }
                                   }
                               } 
                              ch_layout_fix = ch_layout_fix->pointer;
                           }
                       }
                       aux = aux->pointer;
                   }
            }
            if (((errors_summary_old - errors_summary) < (errors_summary * MINIMAL_IMPROVEMENT)) || (errors_summary == 0)) {
                printf("generation %ld: too small progress -> stop\n", generation);
                return 0;
            }
            errors_summary_old = errors_summary;
            coef = coef * 2;

            if ((STEP_SIZE/coef)/2 < 1) {
                printf("generation %ld: The next step is too small %d -> stop\n", generation, (STEP_SIZE/coef)/2);
                return 0;
            } 
        }

        if (generation >= max_generations) {
            printf("generation %ld: it is max_generation -> stop\n", generation);
            return 0;
        }
        generation = generation + 1;
    }

    return 0;
}


int main(int argc, char** argv)
{
    StrList *ligandsnames;
    XYLayout *ligands_layout, *r_ligands_layout;
    XYLayout *sses_layout;
    ThreeDLayout *channels_3d_layout;

    StrList *all_s, *prime_s = NULL;
    AngleLayout *sses_sizes;
    TabDiff *lsdistances;

    int STEP_SIZE;
    int CHECK_EVERY;

    // file openning
    char *filename;
    FILE *fp;

    // file lines reading
    int MaxLineLen = 128;
    char ch;
    char *LineBuffer = (char *)malloc(sizeof(char) * MaxLineLen);
    int LinePos = 0;

    // variable assigning
    int res;
    int counter = 0;

    // cJSON variables
    cJSON *json;               // whole json file
    cJSON *cjson_buffer_var;   // level1 json variable
    cJSON *child;              // auxiliary variable
    int len;
    char *buffer;              // json file buffer

    // file openning
    filename = argv[1];
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }
    fseek(fp, 0L, SEEK_END);
    int length = ftell(fp);
    buffer = malloc((length+1)*sizeof(char));
    fclose(fp);
    
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }
    len = fread(buffer, 1, length+1, fp);
    fclose(fp);

    // parse the JSON data
    json = cJSON_Parse(buffer);
    if (json == NULL) {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr != NULL) {
            printf("Error: %s\n", error_ptr);
        }
        cJSON_Delete(json);
        return 1;
    }

    // read STEP_SIZE variable (int variable)
    // EXAMPLE: "STEP_SIZE":	56
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "STEP_SIZE");
    STEP_SIZE = cjson_buffer_var->valueint;
    //printf("Name: STEP_SIZE -> %d \n", STEP_SIZE);

    // read CHECK_EVERY variable (int variable)
    // EXAMPLE: "CHECK_EVERY":	13900
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "CHECK_EVERY");
    CHECK_EVERY = cjson_buffer_var->valueint;
    //printf("Name: CHECK_EVERY -> %d \n", CHECK_EVERY);

    // read ligandsnames variable (int variable)
    // EXAMPLE: "ligandsnames":	["20.06A0", "20.06A1", "20.06A2", "20.06A3", ....
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "ligandsnames");
    ligandsnames = read_str_list_variable(cjson_buffer_var);
    //print_list("ligandsnames", ligandsnames);

    // read ligandsnames variable (int variable)
    // EXAMPLE: "ligandsnames":	["20.06A0", "20.06A1", "20.06A2", "20.06A3", ....
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "lsdistances");
    lsdistances = read_tabdiff_variable(cjson_buffer_var);
    //print_tabdiff("lsdistances", lsdistances);

    // read ligands_layout variable (int variable)
    // EXAMPLE: "ligands_layout": {"20.06A0":[0, 0],"20.06A1":[0, 0],...
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "ligands_layout");
    ligands_layout = read_xylayout_variable(cjson_buffer_var);
    //print_xylayout("ligands_layout", ligands_layout);

    // read all_s variable 
    // EXAMPLE: "all_s": ['H01AALL', 'H02AALL',  ...
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "all_s");
    all_s = read_str_list_variable(cjson_buffer_var);
    //print_list("all_s", all_s);

    // read prime_s variable 
    // EXAMPLE: "prime_s": ['H01AALL', 'H02AALL',  ...
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "prime_s");
    prime_s = read_str_list_variable(cjson_buffer_var);
    //print_list("prime_s", prime_s);

    // read sses_sizes variable 
    // EXAMPLE: "sses_sizes": ['H01AALL', 'H02AALL',  ...
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "sses_sizes");
    sses_sizes = read_anglelayout_variable(cjson_buffer_var);
    //print_anglelayout("sses_sizes", sses_sizes);

    // read sses_layouts variable
    // EXAMPLE: "sses_layout": 
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "sses_layout");
    sses_layout = read_xylayout_variable(cjson_buffer_var);
    //print_xylayout("sses_layout", sses_layout);

    printf("all data red\n ");

    //printf("We modify ligands positions\n");
    //print_list("ligandsname", ligandsnames);

    res = iteration(ligandsnames, ligands_layout, all_s, prime_s, sses_sizes, sses_layout, lsdistances, STEP_SIZE, CHECK_EVERY);

    if (res != 0) {
        printf("error\n");
    } else {
        fp = fopen(argv[2], "w");
        if (fp == NULL)
        {
            printf("Error: could not open file %s\n", filename);
            return 1;
        }
        r_ligands_layout = ligands_layout;
        while (r_ligands_layout != NULL) {
            fprintf(fp, "layout,%s,%.2f,%.2f\n", r_ligands_layout->name, r_ligands_layout->x, r_ligands_layout->y);
            r_ligands_layout = r_ligands_layout->pointer;
        }
        fclose(fp);

    }
//    oo= time(0);
    printf("konec\n");

}

