#include <stdio.h>
#include <gmp.h>

int main() {
    // Initialize the arbitrary-precision integer
    mpz_t big_num;
    mpz_init(big_num);
    
    // Set the large number from a string
    mpz_set_str(big_num, "3744843080529615909019181510330554205500926021947", 10); // 10 is the base
    
    // Print the number
    gmp_printf("The big number is: %Zd\n", big_num);
    
    // Clear the memory occupied by big_num
    mpz_clear(big_num);
    
    return 0;
}
