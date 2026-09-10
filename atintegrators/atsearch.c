/* Implementation of search algorithms for AT. */


int binarySearch(double *array,double value,int upper,int lower,int nStep){
    int pivot = (int)(lower+upper)/2;
    if ((upper-lower)<=1){
        return lower;
    };
    if (value < array[pivot]){
        upper = pivot;
        nStep+=1;
        return binarySearch(array,value,upper,lower,nStep);
    } else if (value > array[pivot]){
        lower = pivot;
        nStep+=1;
        return binarySearch(array,value,upper,lower,nStep);
    }else{
        return pivot;
    };
};

