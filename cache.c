#include <stdlib.h>
#include <string.h>
#include "cache.h"
#include "kernel.h"
#include "globals.h"

// cache anfang

int nr_of_cache_entries;
int nr_of_elements;

double** data;

int* look_up_table; // translates data id to cache position
int* reverse_look_up_table; // translates cache positions to id
int* circular_array; // safes order in which cache pos was inserted

int ca_first;
int ca_last;
int ca_free_pos; // safes which pos has no yet been occupied

void init(int noce, int noe) {
    nr_of_cache_entries = noce;
    nr_of_elements = noe;

    // allocate memory
    look_up_table = malloc( sizeof(int) * nr_of_elements );
    memset( look_up_table, -1, sizeof(int) * nr_of_elements );

    reverse_look_up_table = malloc( sizeof(int) * nr_of_cache_entries );
    memset( reverse_look_up_table, -1, sizeof(int) * nr_of_cache_entries );

    circular_array = malloc( sizeof(int) * nr_of_cache_entries );
    memset( circular_array, -1, sizeof(int) * nr_of_cache_entries );

    data = malloc( sizeof( double* ) * nr_of_cache_entries);

    int i;
    for(i=0;i<nr_of_cache_entries; i++)
    {
        data[i] = malloc( sizeof( double) * (prob[0].l + prob[1].l) );
    }

    // init pointer
    ca_first = 0;
    ca_last = nr_of_cache_entries - 1;

}

void cleanup() {
    free(reverse_look_up_table);
    free(look_up_table);
    free(circular_array);
    int i;
    for(i=0;i<nr_of_cache_entries; i++)
    {
        free( data[i] );
    }
    free(data);
}

void get_data(int id, int set, int pointer)
{
//    data[pointer] = (double) id * id;

    int i;
    for(i=0;i<prob[0].l;i++)
    {
        //printf("set1 = %d, id = %d,  set2 = %d, id = %d res = %f\n", set, id, 0, i,  kernel(set, id, 0, i));
        data[pointer][i] = kernel(set, id, 0, i);
    }

    for(i=0;i<prob[1].l;i++)
    {
        //printf("set1 = %d, id = %d,  set2 = %d, id = %d   res = %f \n", set, id, 1, i,  kernel(set, id, 1, i));
        data[pointer][i + prob[0].l] = kernel(set, id, 1, i);
    }
}


void ca_add(int id) {
    int last_id = reverse_look_up_table[ circular_array[ca_last] ]; // clean up look up table
    if(circular_array[ca_last] != -1)
    {
        //pos = look_up_table[ last_id ];
        look_up_table[ last_id ] = -1;
    } else {
        circular_array[ca_last] = ca_free_pos;
        ca_free_pos++;
    }

    //circular_array[ca_last] = pos;
    ca_first = ca_last;
    ca_last = ca_last - 1;
    if(ca_last<0) ca_last = nr_of_cache_entries - 1;

    reverse_look_up_table[circular_array[ca_first]] = id;
    look_up_table[id] = circular_array[ca_first];
}

void ca_bring_forward(int pos)
{
//    printf("bring_fordward. enter. pos = %d\n", pos);
    int current = ca_first;
    int pos_temp = circular_array[current];
    int pos_temp2 = -1;
//    int i;
//    printf("circular array: ");
//    for(i=0; i< nr_of_cache_entries; i++)
//        printf(" %d: %d - ", i, circular_array[i]);
//    printf("\n");

//    printf("lut: ");
//    for(i=0; i< nr_of_elements; i++)
//        printf(" %d: %d - ", i, look_up_table[i]);
//    printf("\n");

//    printf("first = %d   last = %d \n", ca_first, ca_last);


    do{
//        printf("bring_fordward. cycle. \n");

        pos_temp2 = pos_temp;
        current = current + 1;
        if(current>=nr_of_cache_entries) current = 0;
        pos_temp = circular_array[current];
//        printf("current = %d   last = %d  pt = %d, pt2 = %d\n", current, last, pos_temp, pos_temp2);
        circular_array[current] = pos_temp2;

//    printf("circular array 2: ");
//    for(i=0; i< nr_of_cache_entries; i++)
//        printf(" %d: %d - ", i, circular_array[i]);
//    printf("\n");

    } while( pos_temp != pos);

    circular_array[ca_first] = pos;

//    printf("circular array 3: ");
//    for(i=0; i< nr_of_cache_entries; i++)
//        printf(" %d: %d - ", i, circular_array[i]);
//    printf("\n");

}

double* get_element(int id, int set)
{
    //printf(" get_element(): id = %d, set = %d \n", id, set);
    int idset = id + set* prob[0].l;
    if( look_up_table[idset] == -1 ) { // cache miss
        ca_add(idset);
        get_data(id, set, circular_array[ca_first]);
        //printf("cache miss, id = %d, set = %d\n", id, set);
    } else { //cache hit
        //printf("cache hit\n");
        if(look_up_table[idset] != circular_array[ca_first])
        {
            ca_bring_forward(look_up_table[idset]);
        }
    }
    //printf("get_element = data[%d]  ca_first = %d \n", circular_array[ca_first], ca_first);
    return data[circular_array[ca_first]];
}

// cache ende

