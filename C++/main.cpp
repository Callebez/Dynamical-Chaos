#include<iostream>
#include<iterator> // for iterators
#include<vector> // for vectors
using namespace std;
int main()
{
    vector<double> ar = { 1.564, 2.54, 3, 4.65464, 5.654 };
      
    // Declaring iterator to a vector
    vector<double>::iterator ptr;
      
    // Displaying vector elements using begin() and end()
    cout << "The vector elements are : ";
    for (ptr = ar.begin(); ptr < ar.end(); ptr++)
        cout << *ptr << " ";
      
    return 0;    
}