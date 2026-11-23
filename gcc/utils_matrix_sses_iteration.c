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
   float angle;
   void *pointer;
} FullLayout;

typedef struct {
   char *name;
   char *cname;
   float x;
   float y;
   void *pointer;
} ClusterLayout;

typedef struct {
   char *name;
   char *cname;
   float angle;
   void *pointer;
} AngleClusterLayout;


typedef struct {
   char *name;
   char *cname;
   float x;
   float y;
   float angle;
   void *pointer;
} FullClusterLayout;


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
} TDLevelOneRow;



typedef struct {
   char *name;
   float row_size;
   void *row;
   TDLevelOneRow *pointer;
} TabDiff;

/*
char *flags;
int STEP_SIZE;
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


char print_fulllayout(char *name, FullLayout* values)
{
    FullLayout *rec = values;
    printf("\n-- %s --: ", name);
    while (rec != NULL) {
        printf("%s<x:%4.1f,y:%4.1f a:%5.1f>, ", rec->name, rec->x, rec->y, rec->angle);
        rec = rec->pointer;
    }
    printf("\n");
}


char print_clusterlayout(char *name, ClusterLayout* values)
{
    ClusterLayout *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s,[%s]<x:%.2f,y:%.2f>, ", rec->name, rec->cname, rec->x, rec->y);
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


char print_fullclusterlayout(char *name, FullClusterLayout* values)
{
    FullClusterLayout *rec = values;
    printf("-- %s --: ", name);
    while (rec != NULL) {
        printf("%s,[%s]<x:%.2f,y:%.2f><angle:%.2f>, ", rec->name, rec->cname, rec->x, rec->y, rec->angle);
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
    TDLevelOneRow *chdiff;
    TabDiff *row;
    printf("-- %s --: \n", name);
    row = values;
    while (row != NULL) {
        chdiff = row->pointer;
        printf("   %4s %p ::\n", row->name, row->pointer);
        while (chdiff != NULL) {
            diff = chdiff->pointer;
            while (diff != NULL) {
                printf("      {%2.0f}x{%2.0f}(%2.0f)(%s) \n", diff->diff_value, diff->twod_diff, diff->size, diff->name);
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
    TDLevelOneRow *chdiff;
    TabDiff *row;
    printf("-- %s --: \n", name);
    row = values;
    while (row != NULL) {
        chdiff = row->pointer;
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




// read input structures
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


ClusterLayout *read_clusterlayout_variable(cJSON *cjson_buffer_var)
{
    ClusterLayout *list = NULL;
    ClusterLayout *rec = NULL, *new_rec;

    cJSON *cluster, *sse, *coords;
    cluster = cjson_buffer_var->child;
    int val;

    while (cluster != NULL) {
        // go through each cluster
        sse = cluster->child;
        while (sse != NULL) {
            // go through each sse in cluster
            coords = sse->child;

            new_rec = (ClusterLayout *)malloc(sizeof(ClusterLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = sse->string;
            new_rec->cname = cluster->string;
            new_rec->x = coords->valuedouble;
            coords = coords->next;
            new_rec->y = coords->valuedouble;
            new_rec->pointer = NULL;
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
            } else {                                        // not the first, connect it to the line
                rec->pointer = new_rec;
            }
            rec = new_rec;
            sse = sse->next;
        }
        cluster = cluster->next;
    }
    return list;
}


AngleClusterLayout *read_angleclusterlayout_variable(cJSON *cjson_buffer_var)
{
    AngleClusterLayout *list = NULL;
    AngleClusterLayout *rec = NULL, *new_rec;

    cJSON *cluster, *sse;
    cluster = cjson_buffer_var->child;
    int val;

    while (cluster != NULL) {
        // go through each cluster
        sse = cluster->child;
        while (sse != NULL) {
            // go through each sse in cluster
            new_rec = (AngleClusterLayout *)malloc(sizeof(AngleClusterLayout)); // create new record rec->name is the name, pointer to null
            new_rec->name = sse->string;
            new_rec->cname = cluster->string;
            new_rec->angle = (float)(sse->valueint);
            new_rec->pointer = NULL;
            if (list == NULL) {                             // the first record, mark it as start
                list = new_rec;
             } else {                                        // not the first, connect it to the line
                rec->pointer = new_rec;
            }
            rec = new_rec;
            sse = sse->next;
        }
        cluster = cluster->next;
    }
    return list;
}


char *sses_channel_prefix(char *name)
{
    int counter = 0;
    while (counter<strlen(name)) {
        if (isalpha(name[counter])) {
            if ((counter+1<strlen(name)) && (isdigit(name[strlen(name)-1]))) {            
                return getSubString(0, counter, name);
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
    TDLevelOneRow *row_levelone_new, *row_levelone_r;

    cJSON *cjson_row, *child;
    cjson_row = cjson_buffer_var->child;
    int res;
    int leave;
    char *prefix;
    char *ret;
    int find;

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
        while (child != NULL) {
            prefix = &(child->string[strlen(child->string)-1]);
            row_levelone_r = row->pointer;    // the first row_levelone record
            // find whether there is level one structure for our diff record
            find = 1;
            while ((find == 1) && (row_levelone_r != NULL)) {
                if (strcmp(row_levelone_r->name, prefix) == 0) {
                    // we find proper structure yet
                    find = 0;
                } else {
                    // we do not find proper structure yet
                    row_levelone_r = row_levelone_r->row;
                }
            }
            if ((row_levelone_r == NULL) || (find == 1)) {
                //there was no level one record suitable for the diff
                // we have to create new one new channel line
                row_levelone_new = (TDLevelOneRow *)malloc(sizeof(TDLevelOneRow));
                row_levelone_new->name = prefix;
                row_levelone_new->pointer = NULL;
                row_levelone_new->row = NULL;

                if (row->pointer == NULL) {
                    // connect tohe row_level one to the row structure (it was the first one)
                    row->pointer = row_levelone_new;
                } else {
                    // connect the row_levelone structure to the previous one structures (there already were some)
                    row_levelone_r = row->pointer; 
                    while (row_levelone_r->row != NULL) {
                        row_levelone_r = row_levelone_r->row;
                    }
                    row_levelone_r->row = row_levelone_new;
                }
                row_levelone_r = row_levelone_new;
            } else {
            }
            diff_new = (Diff *)malloc(sizeof(Diff));  // create new diff instance (name,value, null, null)
            diff_new->name = child->string;
            diff_new->diff_value = child->valueint;
            diff_new->pointer = NULL;
            if (row_levelone_r->pointer == NULL) {
                // start of the new row for new sse
                row_levelone_r->pointer = diff_new;
            } else {
                // we have to find proper position
                diff = row_levelone_r->pointer;
                while (diff->pointer != NULL) {
                    diff = diff->pointer;
                }
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
        coords = coords->next;
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

// =====================================================================================
// all data read
// ====================================================================================


FullLayout *agregate_xy_angle_layout_together(XYLayout *layout, AngleLayout *angle_layout)
{
    FullLayout *list = NULL, *rec = NULL, *rec_new;
    XYLayout *xy_rec=layout;
    AngleLayout *angle_rec;
    while (xy_rec) { 
        rec_new = (FullLayout *)malloc(sizeof(FullLayout));      // create new diff instance (name,value, null, full)
        rec_new->name = (char *)malloc(sizeof(char)*(strlen(xy_rec->name)+1));
        strcpy(rec_new->name, xy_rec->name);
        rec_new->x = xy_rec->x;
        rec_new->y = xy_rec->y;
        rec_new->pointer = NULL;

        angle_rec = angle_layout;                                // angle copy from a_layout
        while (angle_rec) {
            if (strcmp(angle_rec->name, rec_new->name) == 0) {
                 rec_new->angle = angle_rec->angle;
                 break;
            }
            angle_rec = angle_rec->pointer;
        }

        if (list == NULL) {
            list = rec_new;
        } else {
            rec->pointer = rec_new;                              // row only continue
                    }
        rec = rec_new;
        xy_rec = xy_rec->pointer;
    }
    return list;
}

FullClusterLayout *agregate_xy_angle_cluster_layout_together(ClusterLayout *cluster_layout, AngleClusterLayout *angle_cluster_layout)
{
    FullClusterLayout *list = NULL, *rec = NULL, *rec_new = NULL;
    ClusterLayout *c_rec=cluster_layout;
    AngleClusterLayout *angle_rec;

    while (c_rec) {  
        rec_new = (FullClusterLayout *)malloc(sizeof(FullClusterLayout));      // create new diff instance (name,value, null, null)
        rec_new->name = (char *)malloc(sizeof(char) * (strlen(c_rec->name)+1));
        strcpy(rec_new->name,c_rec->name);
        rec_new->cname = (char *)malloc(sizeof(char) * (strlen(c_rec->cname)+1));
        strcpy(rec_new->cname,c_rec->cname);
        rec_new->x = c_rec->x;
        rec_new->y = c_rec->y;
        rec_new->pointer = NULL;
        angle_rec = angle_cluster_layout;                                // angle copy from a_layout
        while (angle_rec) {
            if (strcmp(angle_rec->cname, rec_new->cname) == 0) {
                 rec_new->angle = angle_rec->angle;
                 break;
            }
            angle_rec = angle_rec->pointer;
        }
        if (list == NULL) {
            list = rec_new;
        } else {
            rec->pointer = rec_new;                              // row only continue
                    }
        rec = rec_new;
        c_rec = c_rec->pointer;
    }
    return list;
}


int FindSubstring(char *txt, char *sub) {
    int n = strlen(txt);
    int m = strlen(sub);
    int i, j;

    // Iterate through txt
    for (i = 0; i <= n - m; i++) {
        // Check for substring match
        for (j = 0; j < m; j++) {
            // Mismatch found
            if (txt[i + j] != sub[j]) {
                break;
            }
        }
        // If we completed the inner loop, we found a match
        if (j == m) {
            return i;
        }
    }
    // No match found
    return -1;
}


FullClusterLayout *count_sse_coords_based_on_cluster_coords(FullLayout *full_layout, FullClusterLayout *full_cluster_layout)
{
    FullClusterLayout *list = NULL, *rec = NULL, *rec_new;
    FullClusterLayout *c_res = NULL;
    float angle;
    char *name, *aux;
    int start, counter;
    float x,y;
    float ssedif, beta;
    float coord_sse1_x, coord_sse1_y;
    int l;

    // go through all clusters
    while (full_layout != NULL) {
        angle = full_layout->angle/180 * M_PI;
        name = (char *)malloc(sizeof(char)*(strlen(full_layout->name)+1));
        strcpy(name, full_layout->name);
        l = strlen(name);
        start = 0;
        counter = 0;
        while (counter <= l){
            // we have to go through all sses from cluster
            if ((name[counter] == ':') || (name[counter] == '\0')) {    // we have sse name
               // we found one sse from cluster (its name will be in aux)
               aux = getSubString(start, counter, name);
               start = counter + 1;

               c_res = full_cluster_layout;
               while (c_res != NULL) {
                    // go through all sses in full_cluste_layout to find matching sse record
                    if (strcmp(c_res->name,aux) == 0) {
                        // we find relevant full cluster record for sse aux
                        x = c_res->x;
                        y = c_res->y;
                        ssedif = sqrt(x*x+y*y);
                        if (y == 0) {
                            beta = 1;
                        } else {
                            beta = atan(x/y);
                        }
                        if (y<0) {
                            coord_sse1_x = full_layout->x - (ssedif*sin(angle+beta));
                            coord_sse1_y = full_layout->y - (ssedif*cos(angle+beta));
                        } else {
                            coord_sse1_x = full_layout->x + (ssedif*sin(angle+beta));
                            coord_sse1_y = full_layout->y + (ssedif*cos(angle+beta));
                        }
                        rec_new = (FullClusterLayout *)malloc(sizeof(FullClusterLayout));      // create new diff instance (name,value, null, null)
                        rec_new->name = (char *)malloc(sizeof(char) * (strlen(c_res->name)+1));

                        strcpy(rec_new->name,c_res->name);
                        rec_new->cname = (char *)malloc(sizeof(char) * (strlen(c_res->cname)+1));

                        strcpy(rec_new->cname,c_res->cname);
                        rec_new->x = coord_sse1_x;
                        rec_new->y = coord_sse1_y;
                        rec_new->pointer = NULL;
                        rec_new->angle = ((int)(full_layout->angle + c_res->angle)) % 360;

                        if (list == NULL) {
                            list = rec_new;
                        } else {
                            rec->pointer = rec_new;                              // row only continue
                        }
                        c_res = NULL;
                        rec = rec_new;
                    } else {
                        c_res = c_res->pointer;
                    }
               }
            }
            counter = counter + 1;
        }
        aux = getSubString(start, counter, name);
        full_layout = full_layout->pointer;
    }
    return list;
}


int find_sse_coords(FullClusterLayout *sse_cluster_layout, char *for_helic, float *x, float *y)
{
    FullClusterLayout *rec = sse_cluster_layout;
    int i;
    while (rec) {
        i = FindSubstring(rec ->name, for_helic);
        if ((i >= 0) && (strlen(rec->name) == strlen(for_helic))){
            *x = rec->x;
            *y = rec->y;
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


int get_distance_matrix_2d(FullClusterLayout *sse_cluster_layout, StrList *for_helic_list, TabDiff *tabdiff)
{
    int ret = 0;
    char *name, *aux;
    TabDiff *row = NULL;
    TDLevelOneRow *chdiff;
    Diff *diff = NULL;
    float x1, y1, x2, y2;
    int ret_val;
    StrList *list_member = for_helic_list;


    int start = 0;
    int counter = 0;

    while (list_member != NULL) {
         aux = list_member ->name;
         ret_val = find_sse_coords(sse_cluster_layout, aux, &x1, &y1);
         if (ret_val == 1) {
             printf("Error get_distance_matrix_2d 1 (return value %d)", ret_val);
             return -1;
         }
         start = counter + 1;
         row = tabdiff;
         while (row != NULL) {
            if (strcmp(row->name, aux) == 0) {
              chdiff = row->pointer;
              while (chdiff) {
                diff = chdiff->pointer;
                while (diff) {
                    ret_val = find_sse_coords(sse_cluster_layout, diff->name, &x2, &y2);
                    if (ret_val == 1) {
                        printf("Error get_distance_matrix_2d (return value: %d) name %s \n", ret_val, diff->name);
                    } else {
                        diff->twod_diff = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
                    }
                    diff = diff->pointer;
                }
                chdiff = chdiff->row;
              }
            }
            row = row->row;
         }
         list_member = list_member ->pointer;
    }
    return 0 ;
}


float get_measured_error_2d(TabDiff *tabdiff, StrList *for_helices_list)
{
    TabDiff *rec;
    TDLevelOneRow *chdiff;
    Diff *diff;
    float err = 0, aux;
    StrList *str = for_helices_list;
    int l;

    while (str != NULL) {
      // for all sses in for_helices list  - we have to compute error and sum it up
      l = strlen(str->name);
      rec = tabdiff;
      while (rec) {
        // we have to chose proper row of the table tabdiff
        if (strcmp(rec->name, str->name) == 0) {  
            // this record is for the relevant for_helic_list SSE structure, so go and summ up all records
            chdiff = rec->pointer;
            while (chdiff != NULL) {
                diff = chdiff->pointer;
                while (diff) {
                    // TODO this test can be done during tabdiff generation? or sooner
                    aux = abs(diff->diff_value - diff->twod_diff)*(diff->size+rec->row_size);
                    err = err + aux;
                    diff = diff->pointer; 

                }
                chdiff = chdiff->row;
            }
            rec = NULL;

        } else {
            rec = rec ->row;
        }
      }
      str = str->pointer;
    }
    return(err);
}


float move_error(FullLayout *full_layout, XYLayout *start_layout, StrList *for_helices_list,  AngleLayout *sses_sizes, FullClusterLayout *sse_cluster_layout)
{
    XYLayout *rec_xy = start_layout;
    FullClusterLayout *rec_scl = sse_cluster_layout;
    AngleLayout *rec_al = NULL;
    float err = 0, diff, size;
    while (rec_scl) {
        if (is_string_in_strlist(rec_scl->cname, for_helices_list) == 0) {  // it is in for_helices_list
                size = 0;

                rec_al = sses_sizes;
                while (rec_al) {
                    if (strcmp(rec_al->name, rec_scl->cname) == 0) {
                        size = rec_al->angle;
                        rec_al = NULL;
                    } else {
                        rec_al = rec_al->pointer;
                    }
                }

                rec_xy = start_layout;
                while (rec_xy) {
                    if (strcmp(rec_scl->cname, rec_xy->name) == 0){ 
                        diff = sqrt(((rec_scl->x - rec_xy->x)*(rec_scl->x - rec_xy->x))+ ((rec_scl->y - rec_xy->y)*(rec_scl->y - rec_xy->y)));
                        err = err + diff*size;
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



float count_error(StrList *for_helic_list, FullLayout *full_layout, FullClusterLayout *full_cluster_layout, TabDiff *tabdiff, AngleLayout *sses_sizes, XYLayout *start_layout)
{
    AngleLayout *angle_layout;
    FullClusterLayout *sse_cluster_layout;

    int ret = 0;
    float error = 0;

    sse_cluster_layout = count_sse_coords_based_on_cluster_coords(full_layout, full_cluster_layout);

    ret = get_distance_matrix_2d(sse_cluster_layout, for_helic_list, tabdiff);
    if (ret != 0) {
        printf("Error count_error (return value %d)", ret);
        return (-1);
    }

    error = get_measured_error_2d(tabdiff, for_helic_list);
    if (start_layout != NULL) {
        error = error + move_error(full_layout, start_layout, for_helic_list, sses_sizes, sse_cluster_layout);
    }
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
    TDLevelOneRow *row;
    Diff *diff;
    float err = 0;

    while (rec) {
        // it is in for_helices_list
        rec->row_size = find_sses_size(rec->name, sses_sizes);
        row = rec->pointer;
        while (row) {
            diff = row->pointer;
            while (diff) {
                diff->size = find_sses_size(diff->name, sses_sizes);
                diff = diff->pointer;
            }
            row = row->row;
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

int change_xylayout_with_diffs(FullLayout *full_layout, char *name, float dx, float dy)
{
     FullLayout *res = full_layout;
     while (res) {
         if (strcmp(name, res->name) == 0) {
             res->x = res->x +dx;
             res->y = res->y +dy;
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}


int change_anglelayout_with_diffs(FullLayout *full_layout, char *name, float dangle)
{
     FullLayout *res = full_layout;
     while (res) {
         if (strcmp(name, res->name) == 0) {
             res->angle = res->angle +dangle;
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}


int change_layout_reflection_with_diffs(FullClusterLayout *full_layout, char *name)
{
     FullClusterLayout *res = full_layout;
     while (res) {
         if (strcmp(name, res->cname) == 0) {
             res->x = res->x;
             res->y = - res->y;
             res->angle = 360 - res->angle;
             return (0);
         }
         res = res->pointer;
     }
     return (1);
}


int iteration(StrList *movable_clusters, FullLayout *full_layout, FullClusterLayout *full_cluster_layout, StrList *relevant_sses, TabDiff *tabdiff, AngleLayout *sses_sizes, int STEP_SIZE, int CHECK_EVERY, XYLayout *start_layout)
{
    int ANGLES_SIZE = 45;  // #180#45
    float MINIMAL_IMPROVEMENT = 0.02; //#0.01   # if there was not at least (N*100)% improvement in last CHECK_EVERY generations, call it finished layout
    int max_generations = 60000;

    int generation = 0;
    float coef = 1;
    int r, r1, r2;
    float dx, dy, difangle;
    float error_old, error_new, error_checkpoint = 0;
    int res = 0;
    int counter;

    srand(time(NULL));                      // initialize random function
    int aux_random_number = 0;
    int len_movable_clusters;
    char *for_helic;
    StrList *for_helic_list, *fhl_rec;

    len_movable_clusters = strlist_length(movable_clusters);
    if (STEP_SIZE < 1) {
        STEP_SIZE =1;
    }

    print_tabdiff("333 sses_distances", tabdiff);

    while (generation <= max_generations+1) {
        // choose the moved sse
        r = rand() % len_movable_clusters;
        for_helic = strlist_find_nth_value(r, movable_clusters);
        for_helic_list = char_to_list(for_helic);
        error_old = count_error(for_helic_list, full_layout, full_cluster_layout, tabdiff, sses_sizes, start_layout);
        if (error_old < 0) {
            return (-1);
        }

        // choose the new coords
        dx = ((float)(rand() % ((int)(STEP_SIZE/coef*10))))/10 - STEP_SIZE/(2*coef);
        dy = ((float)(rand() % ((int)(STEP_SIZE/coef*10))))/10 - STEP_SIZE/(2*coef);
        change_xylayout_with_diffs(full_layout, for_helic, dx, dy);

        error_new = count_error(for_helic_list, full_layout, full_cluster_layout, tabdiff, sses_sizes, start_layout);

        if (error_new <-0.1) {
            return (-1);
        }

        if (error_new >= error_old) {
            printf("%d: SSES zustava %s o dx:(%f) dy:(%f) z(%f) err pred(%f) po (%f)\n", generation, for_helic, dx, dy, (STEP_SIZE/coef)/2, error_old, error_new);
            change_xylayout_with_diffs(full_layout, for_helic, -dx, -dy);
        } else {
            error_old = error_new;
            printf("%d: SSES zmena %s o dx:(%f) dy:(%f) z(%f) err pred(%f) po (%f)\n", generation, for_helic, dx, dy, (STEP_SIZE/coef)/2, error_old, error_new);
        }
        // if we move strand try to rotate it - three times
        if (for_helic_list->pointer != NULL) {
          counter = 0;
          if (((int)(ANGLES_SIZE/coef*10)) > 0) {
             while (counter < 3) {
               difangle = ((float)(rand() % ((int)(ANGLES_SIZE/coef*10))))/10 - ANGLES_SIZE/(2*coef);
               change_anglelayout_with_diffs(full_layout, for_helic, difangle);
               error_new = count_error(for_helic_list, full_layout, full_cluster_layout, tabdiff, sses_sizes, start_layout);
               if (error_new <0) {
                   return (-1);
               }

               if (error_new >= error_old) {
                   printf("%d: SSES zustava %s\n", generation, for_helic);
                   change_anglelayout_with_diffs(full_layout, for_helic, -difangle);
               } else {
                   error_old = error_new;
                   printf("%d: SSES zmena rotace %s o difangle:(%f) z (%f) err pred(%f) po (%f)\n", generation, for_helic, difangle, (ANGLES_SIZE/coef)/2, error_old, error_new);
               }
               counter = counter + 1;
             }
          }
        }

        // we have to try reflection (if the cluster is nontrivial and <500 generations
        if ((for_helic_list->pointer != NULL) && (generation < 500)){

            fhl_rec = for_helic_list;

            while (fhl_rec != NULL) {
                change_layout_reflection_with_diffs(full_cluster_layout, fhl_rec->name);
                fhl_rec = fhl_rec->pointer;
            }
            difangle = ((float)(rand() % ((int)(360)))) - 180;
            change_anglelayout_with_diffs(full_layout, for_helic, difangle);

            error_new = count_error(for_helic_list, full_layout, full_cluster_layout, tabdiff, sses_sizes, start_layout);
            if (error_new > error_old) {
                   fhl_rec = for_helic_list;
                   while (fhl_rec != NULL) {
                        change_layout_reflection_with_diffs(full_cluster_layout, fhl_rec->name);
                        fhl_rec = fhl_rec->pointer;
                   }
                   change_anglelayout_with_diffs(full_layout, for_helic, -difangle);
            }
        }

        if ((generation % CHECK_EVERY) == 0) {
            // If there was no significant improvement in last 10k iterations, there is no reason to continue
            error_new = count_error(relevant_sses, full_layout, full_cluster_layout, tabdiff, sses_sizes, start_layout);
            if (generation != 0) {
                if ((error_checkpoint - error_new) < (error_new * MINIMAL_IMPROVEMENT)) {
                    return (0);
                }
            }
            error_checkpoint = error_new;
        }

        if ((STEP_SIZE / coef) > 2) {
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
    TDLevelOneRow *orig_row;
    TDLevelOneRow *new_row = NULL, *new_end_row = NULL;
    Diff *orig_diff;
    Diff *new_diff = NULL, *new_end_diff = NULL;

    while (orig_tabdiff) {
        if (is_string_in_strlist(orig_tabdiff->name, relevant_sses) == 0) {
            // row is in relevant_sses - we have to copy them
            new_tabdiff = (TabDiff *)malloc(sizeof(TabDiff));        // create new row instance (name, null, null)
            new_tabdiff->name = orig_tabdiff->name;
            new_tabdiff->row_size = orig_tabdiff->row_size;
            new_tabdiff->row = NULL;
            new_tabdiff->pointer = NULL;
            if (new_start == NULL) {                             // there was no previous - set is at start
               new_start = new_tabdiff;
            } else {
               new_end->row = new_tabdiff;                              // there was previous - connect the new one to old->row
            }
            new_end = new_tabdiff;
            orig_row = orig_tabdiff->pointer;

            while (orig_row) {
              orig_diff = orig_row->pointer;
              new_row = (TDLevelOneRow *)malloc(sizeof(TDLevelOneRow));        // create new row instance (name, null, null)
              new_row->name = orig_row->name;
              new_row->row = NULL;
              new_row->pointer = NULL;
              if (new_end->pointer == NULL) {                             // there was no previous - set is at start
                 new_end->pointer = new_row;
              } else {
                 new_end_row->row = new_row;                              // there was previous - connect the new one to old->row
              }
              new_end_row = new_row;
              orig_diff = orig_row->pointer;

              while (orig_diff) {
                if (is_string_in_strlist(orig_diff->name, relevant_sses) == 0) {
                    // diff structure is in relevant_sses - we have to copy them
                    new_diff = (Diff *)malloc(sizeof(Diff));        // create new row instance (name, null, null)
//                    strcpy(new_diff->name, orig_diff->name);
                    new_diff->name = orig_diff->name;
                    new_diff->diff_value = orig_diff->diff_value;
                    new_diff->twod_diff = orig_diff->twod_diff;
                    new_diff->size = orig_diff->size;
                    new_diff->pointer = NULL;
                    if (new_end_row->pointer  == NULL) {                             // there was no previous - set is at start
                        new_end_row->pointer = new_diff;
                    } else {
                        new_end_diff->pointer = new_diff;                              // there was previous - connect the new one to old->row
                    }
                    new_end_diff = new_diff;
                }
                orig_diff = orig_diff->pointer;
              }
              orig_row = orig_row->row;
            }
            // now we have to copy diff structures connected to pointers
        }
        orig_tabdiff = orig_tabdiff->row;
    }
    return new_start;
}


int main(int argc, char** argv)
{
    char *flags;
    int STEP_SIZE;
    int CHECK_EVERY;
    StrList *relevant_sses;
    StrList *movable_sses;
    StrList *movable_clusters;
    XYLayout *layout;
    AngleLayout *angle_layout;
    AngleLayout *sses_sizes;
    ClusterLayout *cluster_layout;
    AngleClusterLayout *angle_cluster_layout;
    TabDiff *tabdiff;

    FullLayout *full_layout;
    FullClusterLayout *full_cluster_layout;

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
    XYLayout *res_layout, *start_layout = NULL;
    AngleLayout *res_alayout;
    ClusterLayout *res_clayout;
    AngleClusterLayout *res_aclayout;
    TabDiff *res_tabdiff;

    int res;
    FullLayout *rec_fl;
    FullClusterLayout *rec_fcl;
    int len;

    // cJSON variables
    cJSON *json;                // whole json file
    cJSON *cjson_buffer_var;     //level 1 json variable
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
    printf("length %d\n",length );
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

    printf("po\n");
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

    // read movable_clusters variable
    // EXAMPLE: "movable_clusters":	["20.06A0", "20.06A1", "20.06A2", "20.06A3", ....
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "movable_clusters");
    movable_clusters = read_str_list_variable(cjson_buffer_var);
    print_list("movable_clusters", movable_clusters);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "layout");
    layout = read_xylayout_variable(cjson_buffer_var);
    print_xylayout("layout", layout);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "a_layout");
    angle_layout = read_anglelayout_variable(cjson_buffer_var);
    print_anglelayout("angle_layout", angle_layout);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "c_layout");
    cluster_layout = read_clusterlayout_variable(cjson_buffer_var);
    print_clusterlayout("cluster_layout", cluster_layout);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "ca_layout");
    angle_cluster_layout = read_angleclusterlayout_variable(cjson_buffer_var);
    print_angleclusterlayout("angle_cluster_layout", angle_cluster_layout);

    // read sses_sizes variable 
    // EXAMPLE: "sses_sizes": ['H01AALL', 'H02AALL',  ...
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "sses_sizes");
    sses_sizes = read_anglelayout_variable(cjson_buffer_var);
    print_anglelayout("sses_sizes", sses_sizes);

    // read relevant_sses variable
    // EXAMPLE: "relevant_sses":	["20.06A0", "20.06A1", "20.06A2", "20.06A3", ....
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "relevant_sses");
    relevant_sses = read_str_list_variable(cjson_buffer_var);
//    print_list("relevant_sses", relevant_sses);

    // read movable_sses variable
    // EXAMPLE: "movable_sses":	["20.06A0", "20.06A1", "20.06A2", "20.06A3", ....
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "movable_sses");
    movable_sses = read_str_list_variable(cjson_buffer_var);
//    print_list("movable_sses", movable_sses);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "sses_distance");
    tabdiff = read_tabdiff_variable(cjson_buffer_var);
//    print_tabdiff_2("sses_distances", tabdiff);

    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "start_layout");
    start_layout = read_xylayout_variable(cjson_buffer_var);
//    print_xylayout("start_layout", start_layout);

    // read STEP_SIZE variable (int variable)
    // EXAMPLE: "STEP_SIZE":	56
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "STEP_SIZE");
    STEP_SIZE = cjson_buffer_var->valueint;
//    printf("Name: STEP_SIZE -> %d \n", STEP_SIZE);

    // read CHECK_EVERY variable (int variable)
    // EXAMPLE: "CHECK_EVERY":	13900
    cjson_buffer_var = cJSON_GetObjectItemCaseSensitive(json, "CHECK_EVERY");
    CHECK_EVERY = cjson_buffer_var->valueint;
//    printf("Name: CHECK_EVERY -> %d \n", CHECK_EVERY);

//    printf("------------\n");
//    printf("--- read ---\n");
//    printf("------------\n");
//    print_list("movable_clusters", movable_clusters);
//    print_xylayout("layout", layout);
//    print_anglelayout("angle_layout", angle_layout);
//    print_clusterlayout("cluster_layout", cluster_layout);
//    print_angleclusterlayout("angle_cluster_layout", angle_cluster_layout);
//    print_anglelayout("sses_sizes", sses_sizes);
//    print_list("relevant_sses", relevant_sses);
//    print_list("movable_sses", movable_sses);
//    print_tabdiff("sses_distances", tabdiff);
//    print_xylayout("start_layout", start_layout);
//    printf("-- SSTEP_SIZE --: %i\n", STEP_SIZE); 
//    printf("-- CHECK_EVERY --: %i\n", CHECK_EVERY);

    full_layout = agregate_xy_angle_layout_together(layout, angle_layout);
//    print_fulllayout("full_layout", full_layout);
    full_cluster_layout = agregate_xy_angle_cluster_layout_together(cluster_layout, angle_cluster_layout);

//    print_fullclusterlayout("full_cluster_layout", full_cluster_layout);

    // remove non relevant+sses from tabdiff -> not to bother with them in the next parts (iterations :( TODO
//    print_tabdiff_2("000 sses_distances", tabdiff);
//    print_list("relevant_sses", relevant_sses);
    tabdiff = remove_non_relevant_sses(tabdiff, relevant_sses);
//    print_tabdiff_2("111 sses_distances", tabdiff);

//    print_anglelayout("sses_sizes", sses_sizes);
    res = set_size_to_tabdif(tabdiff, sses_sizes);
    if (res != 0 ) {
        printf("6666\n");
        return (res);
    }
//    print_tabdiff_2("222 sses_distances", tabdiff);

    res = iteration(movable_clusters, full_layout, full_cluster_layout, relevant_sses, tabdiff, sses_sizes, STEP_SIZE, CHECK_EVERY, start_layout);
//    printf("res je %d\n", res);

    if (res != 0) {
        printf("chyba\n");
    } else {
        printf("zapisuju %s\n", argv[2]);
        fp = fopen(argv[2], "w");
        if (fp == NULL)
        {
            printf("Error: could not open file %s", filename);
            return 1;
        }
        rec_fl = full_layout;
        while (rec_fl != NULL) {
            fprintf(fp, "full_layout,%s,%.2f,%.2f,%.2f\n", rec_fl->name, rec_fl->x, rec_fl->y, rec_fl->angle);
            rec_fl = rec_fl->pointer;
        }
        rec_fcl = full_cluster_layout;
        while (rec_fcl != NULL) {
            fprintf(fp, "full_cluster_layout,%s,%s,%.2f,%.2f,%.2f\n", rec_fcl->name, rec_fcl->cname, rec_fcl->x, rec_fcl->y, rec_fcl->angle);
            rec_fcl = rec_fcl->pointer;
        }
        fclose(fp);

    }
}