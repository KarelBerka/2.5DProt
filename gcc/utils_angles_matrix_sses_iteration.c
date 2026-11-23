#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef struct {
    char *name;
    void *pointer;
} StrList;

typedef struct {
   char *name;
   float angle;
   void *pointer;
} AngleLayout;

typedef struct {
   char *name;
   char *cname;
   float angle;
   void *pointer;
} AngleClusterLayout;

typedef struct {
   char name[10];
   float diff_value;
   float twod_diff;
   float size;
   void *pointer;
} Diff;


typedef struct {
   char name[10];
   float row_size;
   void *row;
   Diff *pointer;
} TabDiff;

/*
char *flags;
int ANGLES_SIZE;
int CHECK_EVERY;
StrList *relevant_sses;
StrList *movable_sses;
StrList *movable_clusters;
XYLayout *layout;
AngleLayout *angle_layout;
FullLayout *full_layout;
AngleLayout *sses_sizes;
ClusterLayout *cluster_layout;
AngleClusterLayout *angle_cluster_layout;
TabDiff *tabdiff;
*/
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


// draw structures

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


char print_angleclusterlayout(char *name, AngleClusterLayout* values)
{
    AngleClusterLayout *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s,[%s]<angle:%.2f>, ", rec->name, rec->cname, rec->angle);
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
    TabDiff *row;
    printf("-- %s --: \n", name);
    row = values;
    while (row != NULL) {
//      if (strcmp(row->name ,"E2") == 0) {
        diff = row->pointer;
        printf("==%4s== ", row->name);
        while (diff != NULL) {
            printf("(%dx%dx%d)(%s) ", (int)diff->diff_value, (int)diff->twod_diff,(int)diff->size, diff->name);
            diff = diff->pointer;
        }
        printf("konec \n");
        row = row->row;
//      }
    }
}

// read input structures


char *read_string_variable(char *LineBuffer, char *name)
{
    char *pos;
    char *res;
    int sln;

    sln = strlen(name);
    if (strlen(LineBuffer)>sln) {
        if (strncmp(LineBuffer, name, sln) == 0) {
            pos = strchr(LineBuffer, ' ');
            res = (char *)malloc(sizeof(char) * (strlen(pos)+1));
            strcpy(res ,pos);
            return(res);
        }
    }
    return("");
}


int read_int_variable(char *LineBuffer, char *name)
{
    char *pos;
    int res;
    int sln;

    sln = strlen(name);
    if (strlen(LineBuffer)>sln) {
        if (strncmp(LineBuffer, name, sln) == 0) {
            pos = strchr(LineBuffer, ' ');
            res = (strtol(pos, NULL, 10));
            return(res);
        }
    }
    return 0;
}


StrList *read_str_list_variable(char *LineBuffer, char *name)
{
    char *pos;
    int sln = strlen(name);

    StrList *list = NULL;
    StrList *rec = NULL, *new_rec;

    int start = 0, end = 0;

    // filter only records which starts to name value
    if (strlen(LineBuffer)<=sln) {
        return list;
    }
    if (strncmp(LineBuffer, name, sln) != 0) {
        return list;
    }

    // parse the row
    pos = strchr(LineBuffer, ' ');
    for (int i=0; i<strlen(pos); i++) {
        if (pos[i] == '\'') {
            if (start == 0) {                                   // start string record
               start = i+1;
            } else {
                end = i;                                        // end string record
                new_rec = (StrList *)malloc((sizeof(StrList)+1));   // create new record rec->name is the name, pointer to null
                new_rec->name = getSubString(start, end, pos);
                new_rec->pointer = NULL;
                if (list == NULL) {                             // the first record, mark it as start
                    list = new_rec;
                } else {                                        // not the first, connect it to the line
                    rec->pointer = new_rec;
                }
                rec = new_rec;
                start = 0;
                end = 0;
            }
        }
    }
    return list;
}


AngleLayout *read_anglelayout_variable(char *LineBuffer, char *name)
{
    char *pos;
    int sln= strlen(name);

    AngleLayout *list = NULL;
    AngleLayout *rec = NULL, *new_rec;

    int sse_start = 0, sse_end = 0, val_start = 0;

    if (strlen(LineBuffer) <= sln) {
        return list;
    }
    if (strncmp(LineBuffer, name, sln) != 0) {
        return list;
    }

    pos = strchr(LineBuffer, ' ');
    for (int i=0; i<strlen(pos); i++) {
        if (pos[i] == '\'') {
           if (sse_start == 0) {
               sse_start = i+1;
           } else {
               sse_end = i;
               val_start = i+3;
           }
        }
        if ((pos[i] == ')')  && (val_start != 0)) {
            new_rec = (AngleLayout *)malloc(sizeof(AngleLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = getSubString(sse_start, sse_end, pos);
            new_rec->angle = atof(getSubString(val_start, i, pos));
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
            } else {
                rec->pointer = new_rec;
            }
            rec = new_rec;

            sse_start = 0;
            sse_end = 0;
            val_start = 0;
        }
    }
    return list;
}


AngleClusterLayout *read_angleclusterlayout_variable(char *LineBuffer, char *name)
{
    char *pos;
    int sln;

    AngleClusterLayout *list = NULL;
    AngleClusterLayout *rec = NULL, *new_rec;

    int sse1_start = 0, sse1_end = 0;
    int sse2_start = 0, sse2_end = 0;
    int val_start = 0;

    sln = strlen(name);
    if (strlen(LineBuffer) <= sln) {
        return list;
    }
    if (strncmp(LineBuffer, name, sln) != 0) {
        return list;
    }

    pos = strchr(LineBuffer, ' ');
    for (int i=0; i<strlen(pos); i++) {
        if (pos[i] == '\'') {
           if (sse1_start == 0) {
               sse1_start = i+1;
           } else if (sse1_end == 0) {
               sse1_end = i;
           } else if (sse2_start == 0) {
               sse2_start = i+1;
           } else {
               sse2_end = i;
           }
        }
        if (pos[i] == ':') {
            val_start = i+1;
        }
        if (pos[i] == ')') {
            sse2_start = 0;
            sse2_end = 0;
        }
        if ((pos[i] == '}') && (val_start != 0)) {
            new_rec = (AngleClusterLayout *)malloc(sizeof(AngleClusterLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = getSubString(sse1_start, sse1_end, pos);
            new_rec->cname = getSubString(sse2_start, sse2_end, pos);
            new_rec->angle = atof(getSubString(val_start, i, pos));
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
            } else {                                        // not the first, connect it to the line
                rec->pointer = new_rec;
            }
            rec = new_rec;

            val_start = 0;
            sse2_start = 0;
            sse2_end = 0;
            sse1_start = 0;
            sse1_end = 0;
        }
        if ((pos[i] == ',') && (val_start != 0)) {
            new_rec = (AngleClusterLayout *)malloc(sizeof(AngleClusterLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = getSubString(sse1_start, sse1_end, pos);
            new_rec->cname = getSubString(sse2_start, sse2_end, pos);
            new_rec->angle = atof(getSubString(val_start, i, pos));
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
            } else {                                        // not the first, connect it to the line
                rec->pointer = new_rec;
            }
            rec = new_rec;

            val_start = 0;
            sse2_start = 0;
            sse2_end = 0;
        }
    }

    return list;
}


TabDiff *read_tabdiff_variable(char *LineBuffer, char *name)
{
    // read line: sses_distance {'H3': {'H3': 0.0, 'H4': 10.449864918385309, 'H5': 10.327968099413358, 'H6': 23.451171271388553, 'H7': 18.139931031084433, 'H8': 17.493395353675627, 'H10': 26.46378059253817, 'E11': 28.046328324755812, 'E12': 24.21834831279788}, 'H4': {'H3': 10.449864918
    // to row structure - each outerr sses per one row record, inner sses are in diff's linked to row

    char *pos;            // position in input row which we are just reading
    int sln;              // length of the input name

    TabDiff *row_start = NULL, *row = NULL , *row_new = NULL;
    Diff *diff = NULL, *diff_new = NULL;

    int sse1_start = 0, sse1_end = 0;     // index of sse1 name position in LineBuffer
    char *sse1_name;                      // sse1 name
    int sse2_start = 0, sse2_end = 0;     //index of sse1 name position in LineBuffer
    char *sse2_name;                      // sse2 name
    int value_start = 0, value_end = 0;   // index of value name position in LineBuffer
    char *value_name;                     // value name

    sln = strlen(name);
    printf("**** %s\n", LineBuffer);
    if (strlen(LineBuffer)>sln) {
        if (strncmp(LineBuffer, name, sln) == 0) {
            pos = strchr(LineBuffer, ' ');
            for (int i=0; i<strlen(pos); i++) {
                printf("<%c>", pos[i]);
                if (pos[i] == '\'') {
                   // we read sses stuff (variables sse1_start - start1_end - the first one , sse2_start sse2_end the second one
                   if (sse1_start == 0) {
                       printf("start_of_row");
                       sse1_start = i+1;
                   } else if (sse1_end == 0) {
                       printf("end_of_row");
                       sse1_end = i;
                       sse1_name = getSubString(sse1_start, sse1_end, pos); // sse1_name set

                       row_new = (TabDiff *)malloc(sizeof(TabDiff));        // create new row instance (name, null, null)
                       strcpy(row_new->name, sse1_name);
                       row_new->row = NULL;
                       row_new->pointer = NULL;
                       if (row_start == NULL) {                             // there was no previous - set is at start
                           row_start = row_new;
                       } else {
                           row->row = row_new;                              // there was previous - connect the new one to old->row
                       }
                       row = row_new;
                   } else if (sse2_start == 0) {
                       sse2_start = i+1;
                   } else {
                       sse2_end = i;
                       sse2_name = getSubString(sse2_start, sse2_end, pos); // sse2_name set, waiting for value
                       sse2_start = 0;
                       sse2_end = 0;
                   }
                }
                if (pos[i] == ':') {
                    value_start= i+1;                                       // value_start set
                }
                if ((pos[i] == '}') && (sse1_start != 0) ) {                                        // value end - row is ending
                    printf("end_of_rec1");
                    value_end = i;
                    value_name = getSubString(value_start, value_end, pos); // value_name set

                    diff_new = (Diff *)malloc(sizeof(Diff));  // create new diff instance (name,value, null, null)
                    strcpy(diff_new->name,sse2_name);
                    diff_new->diff_value = atof(value_name);
                    diff_new->pointer = NULL;
                    if (row->pointer == NULL) {                             // start of the new row for new sse
                        row->pointer = diff_new;
                    } else {
                        diff->pointer = diff_new;                           // row only continue
                    }
                    diff = diff_new;

                    value_start = 0;
                    value_end = 0;
                    sse2_start = 0;
                    sse2_end = 0;
                    sse1_start = 0;
                    sse1_end = 0;
                }
                if ((pos[i] == ',') && (value_start!= 0)) {
                    printf("end_of_rec2");
                    value_end = i;
                    value_name =getSubString(value_start, value_end, pos); // value_name set

                    diff_new = (Diff *)malloc(sizeof(Diff));               // create new diff instance (name,value, null, null)
                    strcpy(diff_new->name,sse2_name);
                    diff_new->diff_value = atof(value_name);
                    diff_new->pointer = NULL;

                    if (row->pointer == NULL) {                            // start of the new  row for new sse
                        row->pointer = diff_new;
                    } else {
                        diff->pointer = diff_new;                          // row only continue
                    }
                    diff = diff_new;

                    value_start = 0;
                    value_end = 0;
                    sse2_start = 0;
                    sse2_end = 0;
                }
            }

            return row_start;
        }
    }
    return row_start;
}

// =====================================================================================
// all data read
// ====================================================================================

AngleClusterLayout *count_sse_angles_based_on_cluster_angles(AngleLayout *angle_layout, AngleClusterLayout *cluster_angle_layout)
{
    AngleClusterLayout *list = NULL, *rec = NULL, *rec_new;
    AngleClusterLayout *c_res = NULL;
    AngleLayout *al_res = angle_layout;
    float alpha, beta;
    char *name, *aux;
    int start, counter;
    float ssedif;
    float coord_sse1_x, coord_sse1_y;

    while (al_res != NULL) {
        alpha = al_res->angle;

        name = (char *)malloc(sizeof(char)*(strlen(al_res->name)+1));
        strcpy(name, al_res->name);
        start = 0;
        counter = 0;
        while (counter < strlen(name)) {
            if (name[counter] == ':') {    // we have sse name

               aux = getSubString(start, counter, name);
               start = counter + 1;

               c_res = cluster_angle_layout;
               while (c_res != NULL) {
                    if (strcmp(aux,c_res ->cname) == 0) {
                        beta  = c_res->angle;

                        rec_new = (AngleClusterLayout *)malloc(sizeof(AngleClusterLayout));      // create new diff instance (name,value, null, null)
                        rec_new->name = (char *)malloc(sizeof(char) * (strlen(c_res->name)+1));
                        strcpy(rec_new->name,c_res->name);

                        rec_new->cname = (char *)malloc(sizeof(char) * (strlen(c_res->cname)+1));
                        strcpy(rec_new->cname,c_res->cname);

                        rec_new->angle = ((int)(alpha + beta)) % 360;
                        rec_new->pointer = NULL;
                        if (list == NULL) {
                            list = rec_new;
                        } else {
                            rec->pointer = rec_new;                              // row only continue
                        }
                        c_res = c_res->pointer;
                        rec = rec_new;
                    } else {
                        c_res = c_res->pointer;
                    }
               }
            }
            counter = counter + 1;
        }
        aux = getSubString(start, counter, name);
        c_res = cluster_angle_layout;
        while (c_res != NULL) {
            if (strcmp(aux,c_res ->cname) == 0) {
                beta = c_res->angle;

                rec_new = (AngleClusterLayout *)malloc(sizeof(AngleClusterLayout));      // create new diff instance (name,value, null, null)
                rec_new->name = (char *)malloc(sizeof(char) * (strlen(c_res->name)+1));
                strcpy(rec_new->name,c_res->name);

                rec_new->cname = (char *)malloc(sizeof(char) * (strlen(c_res->cname)+1));
                strcpy(rec_new->cname,c_res->cname);

                rec_new->pointer = NULL;
                rec_new->angle = alpha + beta;

                if (list == NULL) {
                    list = rec_new;
                } else {
                    rec->pointer = rec_new;                              // row only continue
                }
                rec = rec_new;
            }
            c_res = c_res->pointer;
        }
        al_res = al_res->pointer;
    }
    return list;
}


int find_sse_angle(AngleClusterLayout *sse_cluster_angle_layout, char *for_helic, float *alpha)
{
    AngleClusterLayout *rec = sse_cluster_angle_layout;
    while (rec) {
        if (strcmp(rec ->cname, for_helic) == 0) {
            *alpha = rec->angle;
            return 0;
        }
        rec = rec->pointer;
    }
    return 1;
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


int get_distance_matrix_2d(AngleClusterLayout *sse_cluster_angle_layout, StrList *for_helic_list, TabDiff *tabdiff)
{
    int ret = 0;
    char *name, *aux;
    TabDiff *row = NULL;
    Diff *diff = NULL;
    float alpha, beta;
    int ret_val;
    StrList *list_member = for_helic_list;


    int start = 0;
    int counter = 0;

    while (list_member != NULL) {
         aux = list_member ->name;
         ret_val = find_sse_angle(sse_cluster_angle_layout, aux, &alpha);
         if (ret_val == 1) {
             printf("Error get_distance_matrix_2d 1 (return value %d)", ret_val);
             return 1;
         }
         start = counter + 1;
         row = tabdiff;
         while (row != NULL) {
            if (strcmp(row->name, aux) == 0) {
                diff = row->pointer;
                while (diff) {
                    ret_val = find_sse_angle(sse_cluster_angle_layout, diff->name, &beta);
                    if (ret_val == 1) {
                        printf("Error get_distance_matrix_2d (return value: %d) name %s \n", ret_val, diff->name);
                    } else {
                        diff->twod_diff = alpha - beta;
                    }
                    diff = diff->pointer;
                }
            }
            row = row->row;
         }
         list_member = list_member ->pointer;
    }
    return(0);
}



float get_measured_error_2d(TabDiff *tabdiff, StrList *for_helices_list)
{
    TabDiff *rec = tabdiff;
    Diff *diff;
    float err = 0, aux;

    while (rec) {
//        printf("rec->name [%s]", rec->name);
        if (is_string_in_strlist(rec->name, for_helices_list) == 0) {  // it is in for_helices_list
//            printf("<<%s>> ", rec->name);

            diff = rec->pointer;
            while (diff) {
                // TODO this test can be done during tabdiff generation? or sooner
//                if (is_string_in_strlist(diff->name, relevant_sses) == 0) {  // it is in for_helices_list
                aux = abs((float)((int)(diff->diff_value - diff->twod_diff) % 360))*(diff->size+rec->row_size);
//                printf("gme2 %f %f %f %f\n", diff->diff_value, diff->twod_diff, diff->size, rec->row_size);
                err = err + aux;
//                }
                diff = diff->pointer;
            }
        }
        rec = rec ->row;
    }
//    printf("err %f", err);
    return(err);
}


float move_error(AngleLayout *start_layout, StrList *for_helices_list, AngleLayout *sses_sizes, AngleClusterLayout *cluster_angle_layout)
{
    AngleLayout *rec_xy = start_layout;
    AngleClusterLayout *rec_scl = cluster_angle_layout;
    AngleLayout *rec_al = NULL;
    float err = 0, diff, size;
    while (rec_scl) {
        
        if (is_string_in_strlist(rec_scl->cname, for_helices_list) == 0) {  // it is in for_helices_list
//                printf("rec_scl->cname %s", rec_scl->cname);
//                printf("is on list\n");
                size = 0;
                rec_al = sses_sizes;
                while (rec_al) {
//                    printf("<%s,%d>", rec_al->name, strcmp(rec_al->name, rec_scl->cname));
                    if (strcmp(rec_al->name, rec_scl->cname) == 0) {
//                        printf("naslo se %s size %f je to stejne s %s\n", rec_al->name, rec_al->angle, rec_scl->cname);
                        size = rec_al->angle;
                        rec_al = NULL;
                    } else {
                        rec_al = rec_al->pointer;
                    }
                }
//                printf("/n");
//                printf("size is %f\n", size);
                rec_xy = start_layout;
                while (rec_xy) {
//                    printf("porovnavam s %s", rec_xy->name);
                    if (strcmp(rec_scl->cname, rec_xy->name) == 0){
                        diff = abs(rec_scl->angle - rec_xy->angle);
                        err = err + diff*size;
//                        printf("\n== naslo se diff to pridavam %s %f (%f * %f) celkem : %f\n", rec_scl->cname, diff*size, diff, size, err );
                        rec_xy = NULL;
                    } else {
                        rec_xy = rec_xy->pointer;
                    }
                }
        }
        rec_scl = rec_scl->pointer;
    }
    return(err);
}


float count_error(StrList *for_helic_list, AngleLayout *angle_layout, AngleClusterLayout *cluster_angle_layout, TabDiff *tabdiff, AngleLayout *sses_sizes, AngleLayout *start_layout)
{
//    AngleLayout *angle_layout;
    AngleClusterLayout *sse_cluster_angle_layout;

    int ret = 0;
    float error = 0;

    sse_cluster_angle_layout = count_sse_angles_based_on_cluster_angles(angle_layout, cluster_angle_layout);

    ret = get_distance_matrix_2d(sse_cluster_angle_layout, for_helic_list, tabdiff);
    if (ret != 0) {
        printf("Error count_error (return value %d)", ret);
        return (-1);
    }
    error = get_measured_error_2d(tabdiff, for_helic_list);

//    print_list("for_helic_list", for_helic_list);    
    float me = move_error(start_layout, for_helic_list, sses_sizes, sse_cluster_angle_layout);
//    printf("chyby -  %f %f\n", error, me);
    error = error + me;

    return error;
}

float find_sses_size(char *name, AngleLayout *sses_sizes)
{
    AngleLayout *rec = sses_sizes;

    while (rec) {
         if (strcmp(rec->name, name) == 0) {
              return(rec->angle);
         }
         rec = rec->pointer;
    }
    return(0);
}


int set_size_to_tabdif(TabDiff *tabdiff, AngleLayout *sses_sizes)
// add size values to difftab table - 
{
    TabDiff *rec = tabdiff;
    Diff *diff;
    float err = 0;

    while (rec) {
        // it is in for_helices_list
        rec->row_size = find_sses_size(rec->name, sses_sizes);
        diff = rec->pointer;
        while (diff) {
            diff->size = find_sses_size(diff->name, sses_sizes);
            diff = diff->pointer;
        }
        rec = rec ->row;
    }
    return(0);
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


int change_anglelayout_with_diffs(AngleClusterLayout *cluster_angle_layout, char *name, float dangle)
{
     AngleClusterLayout *res = cluster_angle_layout;
     while (res) {
         if (strcmp(name, res->name) == 0) {
             res->angle = res->angle +dangle;
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}

int iteration(StrList *movable_clusters, AngleLayout *angle_layout, AngleClusterLayout *cluster_angle_layout, StrList *relevant_sses, TabDiff *tabdiff, AngleLayout *sses_sizes, int CHECK_EVERY, AngleLayout *start_layout)
{

    int ANGLES_SIZE = 45;  // #180#45
    float MINIMAL_IMPROVEMENT = 0.02; //#0.01   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
    int max_generations = 60000;

    int generation = 0;
    float coef = 1;
    int r;    // auxilary random value
    float difangle;    // the angle differenceof angle which will be tested
    float error_old, error_new, error_checkpoint = 0;
    int res = 0;
    printf("start iterace");

    srand(time(NULL));                      // initialize random function
    int len_movable_clusters = strlist_length(movable_clusters);
    char *for_helic;
    StrList *for_helic_list, *fhl_rec;

    if (ANGLES_SIZE < 1) {
        ANGLES_SIZE = 1;
    }

    while (generation <= max_generations+1) {
        // choose the moved sse
//        printf("\n");
        r = rand() % len_movable_clusters;
        for_helic = strlist_find_nth_value(r, movable_clusters);
        for_helic_list = char_to_list(for_helic);
        error_old = count_error(for_helic_list, angle_layout, cluster_angle_layout, tabdiff, sses_sizes, start_layout);
        if (error_old < 0) {
            return (-1);
        }

        difangle = ((float)(rand() % ((int)(ANGLES_SIZE/coef*10))))/10 - ANGLES_SIZE/(2*coef);
        change_anglelayout_with_diffs(cluster_angle_layout, for_helic, difangle);
        error_new = count_error(for_helic_list, angle_layout, cluster_angle_layout, tabdiff, sses_sizes, start_layout);
        if (error_new < 0) {
            return (-1);
        }

        if (error_new >= error_old) {
            change_anglelayout_with_diffs(cluster_angle_layout, for_helic, -difangle);
        } else {
            printf("%d for helic %s %f > %f\n", generation, for_helic, error_old, error_new);
            error_old = error_new;
        }

        if ((generation % CHECK_EVERY) == 0) {
            // If there was no significant improvement in last 10k iterations, there is no reason to continue
            error_new = count_error(relevant_sses, angle_layout, cluster_angle_layout, tabdiff, sses_sizes, start_layout);
            if (generation != 0) {
                printf("for helic %d %s %f ? %f\n", generation, for_helic, error_checkpoint, error_new);

                if ((error_checkpoint - error_new) < (error_new * MINIMAL_IMPROVEMENT)) {
                    return (0);
                }
            }
            error_checkpoint = error_new;
        }

        if ((ANGLES_SIZE / coef) > 2) {
            if (((len_movable_clusters<6) && (generation % 1000 == 0)) || ((len_movable_clusters<12) && (generation % 2000 == 0)) || ((len_movable_clusters<25) && (generation % 4000 == 0)) || ((len_movable_clusters<100) && (generation % 6000 == 0))  || (( generation % 8000 == 0))) {
                coef = coef * 2;
            }
        }

        if (generation >= max_generations) {
            return (0);
        }

        generation = generation + 1;
    }
}


TabDiff *remove_non_relevant_sses(TabDiff *tabdiff, StrList *relevant_sses)
{
    // copy only relevant sses from tabdiff
    TabDiff *orig_tabdiff = tabdiff;
    TabDiff *new_tabdiff = NULL, *new_start = NULL, *new_end = NULL;
    Diff *orig_diff;
    Diff *new_diff = NULL, *new_end_diff = NULL;

    while (orig_tabdiff) {
        if (is_string_in_strlist(orig_tabdiff->name, relevant_sses) == 0) {
            // row is in relevant_sses - we have to copy them
            new_tabdiff = (TabDiff *)malloc(sizeof(TabDiff));        // create new row instance (name, null, null)
            strcpy(new_tabdiff->name, orig_tabdiff->name);
            new_tabdiff->row_size = orig_tabdiff->row_size;
            new_tabdiff->row = NULL;
            new_tabdiff->pointer = NULL;
            if (new_start == NULL) {                             // there was no previous - set is at start
               new_start = new_tabdiff;
            } else {
               new_end->row = new_tabdiff;                              // there was previous - connect the new one to old->row
            }
            new_end = new_tabdiff;
            orig_diff = orig_tabdiff->pointer;
            while (orig_diff) {

                if (is_string_in_strlist(orig_diff->name, relevant_sses) == 0) {
                    // diff structure is in relevant_sses - we have to copy them
                    new_diff = (Diff *)malloc(sizeof(Diff));        // create new row instance (name, null, null)
                    strcpy(new_diff->name, orig_diff->name);
                    new_diff->diff_value = orig_diff->diff_value;
                    new_diff->twod_diff = orig_diff->twod_diff;
                    new_diff->size = orig_diff->size;
                    new_diff->pointer = NULL;
                    if (new_end->pointer  == NULL) {                             // there was no previous - set is at start
                        new_end->pointer = new_diff;
                    } else {
                        new_end_diff->pointer = new_diff;                              // there was previous - connect the new one to old->row
                    }
                    new_end_diff = new_diff;
                }
                orig_diff = orig_diff->pointer;
            }
            // now we have to copy diff structures connected to pointers
        }
        orig_tabdiff = orig_tabdiff->row;

    }
    return new_start;
}


int main(int argc, char** argv)
{
    int CHECK_EVERY; //
    StrList *relevant_sses, *movable_cl; //
    StrList *movable_clusters; //
    AngleLayout *angle_layout, *rec_a; //
    AngleLayout *sses_sizes; //
    AngleLayout *start_layout = NULL; //
    AngleClusterLayout *cluster_angle_layout, *rec_ca; //
    TabDiff *tabdiff;

    // file openning
    char *filename;
    FILE *fp;

    //file lines reading
    int MaxLineLen = 128;
    char ch;
    char *LineBuffer = (char *)malloc(sizeof(char) * MaxLineLen);
    int LinePos = 0;

    //variable assigning
    char *res_chr;
    int res_int;
    StrList *res_p_chr;
    AngleLayout *res_alayout;
    AngleClusterLayout *res_aclayout;
    TabDiff *res_tabdiff;

    int res;

    // file openning
    printf("start 0");
    filename = argv[1];
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }

    printf("start 1");
    // read one character at a time and
    // display it to the output
    while ((ch = fgetc(fp)) != EOF) {
        // read file lines
        if (ch != '\n') {
            if (LinePos >= MaxLineLen-2) {
                char *tmp = realloc( LineBuffer, MaxLineLen *2 *sizeof(char));
                if ( !tmp ) {
                    free( LineBuffer );
                    LineBuffer = NULL;
                    break;
                    printf("Error 5 - ");
                }
                else
                {
                    LineBuffer = tmp;
                    MaxLineLen = MaxLineLen *2;
                }
            }
            LineBuffer[LinePos] = ch;
            LinePos = LinePos + 1;

        } else {
            LineBuffer[LinePos] = '\0';
            LinePos = 0;
            res_int = read_int_variable(LineBuffer, "CHECK_EVERY"); //
            if (res_int != 0) {
                CHECK_EVERY = res_int;
            }
            res_p_chr = read_str_list_variable(LineBuffer, "relevant_sses"); //
            if (res_p_chr != NULL) {
                relevant_sses = res_p_chr;
            }
            res_p_chr = read_str_list_variable(LineBuffer, "movable_cl"); //
            if (res_p_chr != NULL) {
                movable_cl = res_p_chr;
            }
            res_p_chr = read_str_list_variable(LineBuffer, "movable_clusters"); //
            if (res_p_chr != NULL) {
                movable_clusters = res_p_chr;
            }
            res_alayout = read_anglelayout_variable(LineBuffer, "a_layout"); //
            if (res_alayout != NULL) {
                angle_layout = res_alayout;
            }
            res_alayout = read_anglelayout_variable(LineBuffer, "sses_sizes"); //
            if (res_alayout != NULL) {
                sses_sizes = res_alayout;
            }
            res_aclayout = read_angleclusterlayout_variable(LineBuffer, "ca_layout"); //
            if (res_aclayout != NULL) {
                cluster_angle_layout = res_aclayout;
            }
            res_tabdiff = read_tabdiff_variable(LineBuffer, "sses_angles");
            if (res_tabdiff != NULL) {
                tabdiff = res_tabdiff;
            }
            res_alayout = read_anglelayout_variable(LineBuffer, "start_angles"); //
            if (res_alayout != NULL) {
                start_layout = res_alayout;
            }

            //sses_angles 

            //start_angles
        }
    }

    // close the file
    fclose(fp);
    printf("------------\n");
    printf("--- read ---\n");
    printf("------------\n");
//    print_list("movable_clusters", movable_clusters);
//    print_anglelayout("angle_layout", angle_layout);
//    print_angleclusterlayout("angle_cluster_layout", cluster_angle_layout);
//    print_anglelayout("sses_sizes", sses_sizes);
//    return(0);
//    print_list("relevant_sses", relevant_sses);
//    print_list("movable_cl", movable_cl);
//    print_tabdiff("sses_distances", tabdiff);
//    printf("-- CHECK_EVERY --: %i\n", CHECK_EVERY);
//    print_anglelayout("start_layout", start_layout);
//    return(0);

    printf("----------------\n");
    printf("- full layouts -\n");
    printf("----------------\n");

//    print_tabdiff("000 sses_distances", tabdiff);
//    print_list("relevant_sses", relevant_sses);
    tabdiff = remove_non_relevant_sses(tabdiff, relevant_sses);
//    print_tabdiff("111 sses_distances", tabdiff);

    res = set_size_to_tabdif(tabdiff, sses_sizes);
    // TODO test res
//    print_tabdiff("sses_distances", tabdiff);

    printf("-----------\n");
    printf("- iterace -\n");
    printf("-----------\n");

    res = iteration(movable_clusters, angle_layout, cluster_angle_layout, relevant_sses, tabdiff, sses_sizes, CHECK_EVERY, start_layout);

    if (res != 0) {
        printf("chyba\n");
    } else {
        fp = fopen(argv[2], "w");
        if (fp == NULL)
        {
            printf("Error: could not open file %s", filename);
            return 1;
        }
        rec_a = angle_layout;
        while (rec_a != NULL) {
            fprintf(fp, "angle_layout,%s,%.2f,%.2f,%.2f\n", rec_a->name, rec_a->angle);
            rec_a = rec_a->pointer;
        }
        rec_ca = cluster_angle_layout;
        while (rec_ca != NULL) {
            fprintf(fp, "angle_cluster_layout,%s,%s,%.2f,%.2f,%.2f\n", rec_ca->name, rec_ca->cname, rec_ca->angle);
            rec_ca = rec_ca->pointer;
        }
        fclose(fp);

    }

}