#include <iostream>
#include "BigNum.h"

using namespace std;

int main() {
    BigNum num1;
    BigNum num2;

    cin >> num1 >> num2;

    cout << "num1 + num2: " << num1 + num2 << endl;
    cout << "num1 - num2: " << num1 - num2 << endl;
    cout << "num1 * num2: " << num1 * num2 << endl;
    cout << "num1 / num2: " << num1 / num2 << endl;
    cout << "num1+num2*num1-num2: " << num1+num2*num1-num2 << endl;
    cout << "num1+num2+num1*num2: " << num1+num2+num1*num2 << endl;
}