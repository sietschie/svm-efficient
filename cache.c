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
    int last_id = reverse_look_up_table[ ca_last ]; // clean up look up table
    if(last_id != -1)
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

    reverse_look_up_table[ca_first] = id;
    look_up_table[id] = ca_first;
}

void ca_bring_forward(int pos)
{
    int temp = circular_array[pos];
    int pos_last = pos;
    int pos_current = pos - 1;
    if(pos_current < 0) pos_current = nr_of_cache_entries - 1;
    while(pos_last != ca_first)
    {
        circular_array[pos_last] = circular_array[pos_current];

        pos_last = pos_current;
        pos_current--;
        if(pos_current < 0) pos_current = nr_of_cache_entries - 1;
    }

    circular_array[ca_first] = temp;


}

double* get_element(int id, int set)
{
    //printf(" get_element(): id = %d, set = %d \n", id, set);
    int idset = id + set* prob[0].l;
    if( look_up_table[idset] == -1 ) { // cache miss
        ca_add(idset);
        get_data(id, set, circular_array[ca_first]);
//        printf("cache miss, id = %d, set = %d\n", id, set);
    } else { //cache hit
        if(look_up_table[idset] != ca_first)
        {
            ca_bring_forward(look_up_table[idset]);
        }
    }
    return data[circular_array[ca_first]];
}

// cache ende

