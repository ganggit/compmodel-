//# modified by Bob Shine*/

/*Note that this code does not check that left does not exceed
the array bound. You need to add this check, before performing the
swaps - both the one in the loop and the final one outside the loop.*/

#include <iostream>
#define SWAP (a,b) { int t; t=a; a=b; b=t; }
int partition (int num[], int low, int high);
void quicksort(int num[], int low, int high );

using namespace std;

int main ()
{
   int num[11] = {200, 99, 73, 60, 40, 50, 30, 20,10,9,1};
   int i, low, high;
   low = 0; high = 10;
   quicksort (num, low, high);
   for (i=0; i < 11; i++)
      cout << num[i] << endl;
   return 0;
}

int partition(int a[], int low, int high )
   {
   int pivot, left, right, tmp;
   int pivot_item;
   pivot_item = a[low];
   pivot = left = low;
   right = high;
   while ( left < right ) {
     /* Move left while item < pivot */
     while( a[left] <= pivot_item ) left++;
     /* Move right while item > pivot */
     while( a[right] > pivot_item ) right--;
     if ( left < right ) 
         {tmp = a[left]; a[left] = a[right]; a[right] = tmp;}
     }
   /* right is final position for the pivot */
   a[low] = a[right];
   a[right] = pivot_item;
   return right;
   }

void quicksort(int a[], int low, int high )
   {
   int pivot;
   /* Termination condition! */
   if ( high > low )
     {
     pivot = partition( a, low, high );
     quicksort( a, low, pivot-1 );
     quicksort( a, pivot+1, high );
     }
   }
