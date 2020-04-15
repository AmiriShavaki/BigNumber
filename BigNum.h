#ifndef BIGNUM_H
#define BIGNUM_H

#include <iostream>
#include <string>
#include <algorithm> //To use max function

using namespace std;

const int MAXFRACSIZE = 20;

class BigNum {
    //Let us using help of these operator for assigning input big num in our BigNum object
    friend istream& operator>> (istream& input, BigNum& storage);
    friend ostream& operator<< (ostream& output, const BigNum& storage);
    
    friend int getDigit(const BigNum& given); //Give us integer digit of given 1digit BigNum

    private:
    int* integerPart;
    int* fractionalPart;
    bool isPositive; //true -> positive ... false -> negative
    int integerPartSize, fractionalPartSize; //We store sizes because of that later in some operations (e.g print output) they are useful to us

    public:
    BigNum operator= (const BigNum& value);
    BigNum(){}; //Let compiler ignore declearing empty BigNum object
    BigNum(int integerPartSize, int fractionalPartSize, bool isPositive = true);

    BigNum operator+ (const BigNum& value) const;
    BigNum operator- (const BigNum& value) const;
    BigNum operator* (const BigNum& value) const;
    BigNum operator/ (const BigNum& value) const;
    
    //Calculates absolute greater value between two given input
    //Notice that this function doesn't return a > b but return |a| > |b| instead
    bool operator> (const BigNum& value) const; 
    bool operator== (const BigNum& value) const; 
    bool operator>= (const BigNum& value) const; 

    void app(const int digit); //Multiply in 10 then add digit
    const bool isZero() const; //Return true if and only if number equals to zero
    void makeValue(int value); //Usefull to make 0, 1 or 10 value in BigNum objects
    void repairSize(); //Remove junk elements if there is any
};

bool BigNum::operator>= (const BigNum& value) const {
    return (*this > value) || (*this == value);
}

bool BigNum::operator== (const BigNum& value) const {
    return !(*this > value) && !(value > *this);
}

void BigNum::repairSize() {
    int delCnt = 0;
    for (int i = 0; i < integerPartSize - 1; i++) {
        if (integerPart[i] == -1 || integerPart[i] == 0) {
            delCnt++;
        } else {
            break;
        }
    }
    int size = integerPartSize - delCnt;
    int* cpy = new int [size];
    for (int i = 0, j = delCnt; i < size; i++, j++) {
        cpy[i] = integerPart[j];
    }
    delete [] integerPart;
    integerPartSize = size;
    integerPart = cpy;
    fractionalPartSize = min(MAXFRACSIZE, fractionalPartSize);
}

int getDigit(const BigNum& given) {
    return given.integerPart[0];
}

const bool BigNum::isZero() const {
    for (int i = 0; i < integerPartSize; i++) {
        if (integerPart[i] > 0) {
            return false;
        }
    }

    for (int i = 0; i < fractionalPartSize; i++) {
        if (fractionalPart[i] > 0) {
            return false;
        }
    }    

    return true;
}

void BigNum::app(const int digit) {
    if (isZero()) { //Zero case
        integerPart[0] = digit;
    } else { //Non-zero case
        //Copy current value in another dynamic array
        int* tmp = new int [integerPartSize];
        for (int i = 0; i < integerPartSize; i++) {
            tmp[i] = integerPart[i];
        }

        //Replacing new dynamic array
        delete [] integerPart;
        integerPart = new int [integerPartSize + 1];
        for (int i = 0; i < integerPartSize; i++) {
            integerPart[i] = tmp[i];
        }
        integerPart[integerPartSize] = digit;
        integerPartSize++;
    }
}

void BigNum::makeValue(int value) {
    isPositive = true;
    if (value == 10) { //Case of 10
        integerPartSize = 2;
        fractionalPartSize = 0;
        integerPart = new int [2];
        integerPart[0] = 1;
        integerPart[1] = 0;
    } else { //Case of 0 or 1
        integerPartSize = 1;
        fractionalPartSize = 0;
        integerPart = new int [1];
        integerPart[0] = value;
    }
}

BigNum BigNum::operator/ (const BigNum& value) const {
    //If both of input argumants had same signs then result number will be positive, otherwise negative
    bool sign = isPositive == value.isPositive; //false -> negative  true -> positive

    int dividedLength = (this -> fractionalPartSize) + (this -> integerPartSize);
    int divisorLength = value.fractionalPartSize + value.integerPartSize;

    //Delete floating points to make our job easier, later we will add them in satisified positions
    int* divided = new int [dividedLength + MAXFRACSIZE * 2];
    int* res = new int [dividedLength + MAXFRACSIZE * 2]; //As Ms.Baradaran said in Q&A we need to calculate up to 20 extended number in fractional part

    //Copy from *this to divided
    for (int i = 0; i < this -> integerPartSize; i++) {
        divided[i] = (this -> integerPart)[i];
    }
    for (int i = 0; i < this -> fractionalPartSize; i++) {
        divided[i + (this -> integerPartSize)] = (this -> fractionalPart)[i];
    } 
    for (int i = dividedLength; i < dividedLength + MAXFRACSIZE * 2; i++) {
        divided[i] = 0;
    }

    //Copy divisor into a BigNum object
    BigNum divisor(divisorLength, 0);
    for (int i = 0; i < value.integerPartSize; i++) {
        divisor.integerPart[i] = (value.integerPart)[i];
    }
    for (int i = 0; i < value.fractionalPartSize; i++) {
        divisor.integerPart[i + (value.integerPartSize)] = (value.fractionalPart)[i];
    } 
    divisor.repairSize(); //If size didn't compare as well then output of operator> will not be what should be

    //Usefull amounts in next stage
    BigNum zero;
    zero.makeValue(0);
    BigNum ten;
    ten.makeValue(10);
    BigNum one;
    one.makeValue(1);

    //Calculating
    bool flg = false; //Show that if any non-zero digit is placed on answer or not
    BigNum tmp(1, 0);
    int resIndex = 0, floatingPointPos = -(this -> fractionalPartSize) + value.fractionalPartSize + dividedLength;
    for (int i = 0; i < dividedLength + MAXFRACSIZE * 2; i++) { 
        tmp.app(divided[i]);
        while (divisor > tmp) { //Make sure we can divide tmp by divisor (Overloaded operator)
            res[resIndex++] = 0;
            tmp.app(divided[++i]);
        }
        flg = true; //First time we reach here flg changed to true
        int j2 = 0;
        for (BigNum j = zero; ten > j; j = j + one) { //All things here are from BigNumber material!
            if (j * divisor > tmp) { //Overloaded operator
                res[resIndex++] = getDigit(j - one);
                tmp = tmp - ((j - one) * divisor);
                tmp.isPositive = true; //Preventing cause -0 in reminder
                break;
            }

            if (j2 == 9) { // Nine case
                res[resIndex++] = getDigit(j);
                tmp = tmp - (j * divisor);
                tmp.isPositive = true; //Preventing cause -0 in reminder
            }

            j2++; //Update int version of j2
        }
    }

    BigNum ans(floatingPointPos, dividedLength + MAXFRACSIZE * 2, sign);

    //Integer part
    for (int i = 0; i < floatingPointPos; i++) {
        ans.integerPart[i] = res[i];
    } 
    //Fractional part
    for (int i = floatingPointPos, j = 0; i < dividedLength + MAXFRACSIZE * 2; i++, j++) { //j -> index of ans, i -> index of res
        ans.fractionalPart[j] = res[i];
    }

    ans.repairSize();

    return ans;
}

BigNum BigNum::operator* (const BigNum& value) const {
    //Special case: at least one of given BigNumbers equals to zero then we have to return zero
    if (isZero() || value.isZero()) {
        BigNum ans;
        ans.makeValue(0);
        return ans;
    }
    //If both of input argumants had same signs then result number will be positive, otherwise negative
    bool sign = isPositive == value.isPositive; //false -> negative  true -> positive

    int num1Length = (this -> fractionalPartSize) + (this -> integerPartSize), num2Length = value.fractionalPartSize + value.integerPartSize;

    //Delete floating points to make our job easier, later we will add them in satisified positions
    int* num1 = new int [num1Length];
    int* num2 = new int [num2Length];
    int* res = new int [num1Length + num2Length];

    //Copy from *this to num1
    for (int i = 0; i < this -> integerPartSize; i++) {
        num1[i] = (this -> integerPart)[i];
    }
    for (int i = 0; i < this -> fractionalPartSize; i++) {
        num1[i + (this -> integerPartSize)] = (this -> fractionalPart)[i];
    }

    //Copy from value to num2
    for (int i = 0; i < value.integerPartSize; i++) {
        num2[i] = (value.integerPart)[i];
    }
    for (int i = 0; i < value.fractionalPartSize; i++) {
        num2[i + (value.integerPartSize)] = (value.fractionalPart)[i];
    }

    //Initializing res elements to zero
    for (int i = 0; i < num1Length + num2Length; i++) {
        res[i] = -1;
    }

    //Calculating
    for (int i = num2Length - 1; i >= 0; i--) { //multipling num2 digit by digit in num1
        int digit2 = num2[i];
        for (int j = num1Length - 1; j >= 0 ; j--) { //Digit by digit of num1
            int digit1 = num1[j];

            if (res[i + j + 1] == -1) { //For the first time we see this position
                res[i + j + 1] = 0;
            }

            res[i + j + 1] += digit1 * digit2;
        }
    }

    int ansFractionalPartSize = this -> fractionalPartSize + value.fractionalPartSize;
    int ansIntegerPartSize = num1Length + num2Length - ansFractionalPartSize;
    BigNum ans(ansIntegerPartSize, ansFractionalPartSize, sign);

    //Proccessing carries and assigning ans
    for (int i = num1Length + num2Length - 1; i > 0 ; i--) {
        if (res[i - 1] == -1) { //For the first time we see this position
            res[i - 1] = 0;
        }

        res[i - 1] += res[i] / 10;
        res[i] %= 10;
        if (i >= ansIntegerPartSize) {
            ans.fractionalPart[i - ansIntegerPartSize] = res[i];
        } else {
            ans.integerPart[i] = res[i];
        }
    } ans.integerPart[0] = res[0];

    ans.repairSize();

    return ans;
}

//Calculates absolute greater value between two given input
//Notice that this function doesn't return a > b but return |a| > |b| instead
bool BigNum::operator> (const BigNum& value) const {
    if (integerPartSize > value.integerPartSize) { //Check length of integer part
        return true;
    }
    if (integerPartSize < value.integerPartSize) { //Check length of integer part
        return false;
    }
    for (int i = 0; i < integerPartSize; i++) { //Check digit by digit of integer part
        if (integerPart[i] > value.integerPart[i]) { 
            return true;
        }
        if (integerPart[i] < value.integerPart[i]) {
            return false;
        }
    }
    for (int i = 0; i < max(fractionalPartSize, value.fractionalPartSize); i++) { //Check digit by digit of fractional part
        if (fractionalPartSize <= i) { 
            return false;
        }
        if (value.fractionalPartSize <= i) {
            return true;
        }
        if (fractionalPart[i] > value.fractionalPart[i]) { 
            return true;
        }
        if (fractionalPart[i] < value.fractionalPart[i]) {
            return false;
        }        
    }

    return false; //Equal case
}

BigNum BigNum::operator+ (const BigNum& value) const {
    //Maximum size of each part of answer is considered for answer object
    BigNum ans(max(integerPartSize, value.integerPartSize) + 1, max(fractionalPartSize, value.fractionalPartSize), isPositive);
    BigNum cpy = value;

    if (value.isPositive == isPositive) { //If two operands had same signs we want to calculate one of a + b or -(a + b)
        //Fractional part summation
        int carry = 0;
        for (int i = ans.fractionalPartSize - 1; i >= 0; i--) { //i -> ans
            int res = carry; //Result of summation of current digit position plus carry from last position
            if (i < fractionalPartSize) {
                res += fractionalPart[i];
            }
            if (i < value.fractionalPartSize) {
                res += value.fractionalPart[i];
            }
            carry = res / 10;
            ans.fractionalPart[i] = res % 10;
        }

        //Integer part summation
        for (int i = ans.integerPartSize - 1, j = integerPartSize - 1, k = value.integerPartSize - 1; //i -> ans   j -> *this   k -> value
            i >= 0; i--, j--, k--) {
            int res = carry; //Result of summation of current digit position plus carry from last position
            bool flg = res != 0; //Flg shows us this digit is useless to store or not
            if (j >= 0) {
                res += integerPart[j];
                flg = true;
            }
            if (k >= 0) {
                res += value.integerPart[k];
                flg = true;
            }
            
            if (flg) {
                carry = res / 10;
                ans.integerPart[i] = res % 10;
            } else { //If -1 was placed somewhere it means that the digit was considered for carry but currently is useless
                ans.integerPart[i] = -1;
            }
        }
    } else { //If two operands had not same signs we want to calculate one of a - b or b - a by aid of minus operator

        if (isPositive) {
            cpy = value;
            cpy.isPositive = true;
            ans = *this - cpy;
        } else {
            cpy = *this;
            cpy.isPositive = true;
            ans = value - *this;
        }
    }
    ans.repairSize();
    return ans;
}

BigNum BigNum::operator- (const BigNum& value) const {
    bool sign = (*this > value) == isPositive; //false -> negative   true -> positive
    BigNum ans(max(integerPartSize, value.integerPartSize), max(fractionalPartSize, value.fractionalPartSize), sign);

    if (value.isPositive != isPositive) { // a - (-b) or (-a) - b these two cases can be calculated by aid of plus operator
        BigNum cpy = value;
        if (isPositive) { //a - (-b)
            cpy.isPositive = true;
            ans = *this + cpy;
        } else { // (-a) - b
            cpy.isPositive = false;
            ans = *this + cpy;
        }
    } else { //Main part of subtraction operator function
        //Before calculating difference first we should know which one is greater which one is lower
        const BigNum* greater;
        const BigNum* lower; 
        if (*this > value) {
            greater = this;
            lower = &value;
        } else {
            greater = &value;
            lower = this;
        }
        int carry = 0;
        for (int i = max(fractionalPartSize, value.fractionalPartSize) - 1; i >= 0; i--) { //Calculate from the end of fractional part
            if (i >= greater -> fractionalPartSize) {
                ans.fractionalPart[i] = -(lower -> fractionalPart[i]) - carry;
            } else if (i >= lower -> fractionalPartSize) {
                ans.fractionalPart[i] = greater -> fractionalPart[i] - carry;
            } else {
                ans.fractionalPart[i] = greater -> fractionalPart[i] - lower -> fractionalPart[i] - carry;
            }

            if (ans.fractionalPart[i] < 0) {
                ans.fractionalPart[i] += 10;
                carry = 1;
            } else {
                carry = 0;
            }
        }

        //Carry is not deleted from previous loop and we hold it's value to use in next loop

        //Calculate from the end of integer part
        for (int i = (greater -> integerPartSize) - 1, j = (lower -> integerPartSize) - 1; i >= 0; i--, j--) { 
            if (j < 0) {
                ans.integerPart[i] = (greater -> integerPart)[i] - carry;
            } else {
                ans.integerPart[i] = (greater -> integerPart)[i] - (lower -> integerPart)[j] - carry;
            }

            if (ans.integerPart[i] < 0) {
                ans.integerPart[i] += 10;
                carry = 1;
            } else {
                carry = 0;
            }
        }       
    }
    ans.repairSize();
    return ans;
}

BigNum::BigNum(int integerPartSize, int fractionalPartSize, bool isPositive) {
    //Create dynamic arrays fit to data
    this -> integerPartSize = integerPartSize;
    this -> fractionalPartSize = fractionalPartSize;
    this -> isPositive = isPositive;
    integerPart = new int [integerPartSize];
    fractionalPart = new int [fractionalPartSize];
    
    //Initialize all values to zero
    for (int i = 0; i < integerPartSize; i++) {
        integerPart[i] = 0;
    }
    for (int i = 0; i < fractionalPartSize; i++) {
        fractionalPart[i] = 0;
    }
}

BigNum BigNum::operator= (const BigNum& value) {
    //Create new dynamic arrays fit to new data
    integerPart = new int [value.integerPartSize];
    fractionalPart = new int [value.fractionalPartSize];
    
    //Store new sizes
    integerPartSize = value.integerPartSize;
    fractionalPartSize = value.fractionalPartSize;
    
    //Copy process
    for (int i = 0; i < value.integerPartSize; i++) {
        integerPart[i] = value.integerPart[i];
    }
    for (int i = 0; i < value.fractionalPartSize; i++) {
        fractionalPart[i] = value.fractionalPart[i];
    }
    
    isPositive = value.isPositive;

    return value;
}

istream& operator>> (istream& input, BigNum& storage) {
    string raw;
    input >> raw;

    if (raw[0] == '-') {
        storage.isPositive = false;
    } else {
        storage.isPositive = true;
    }

    //Calculate size of integer and fractional part of input number
    storage.integerPartSize = 0;
    //This loop calculates size of integer part
    for (; raw[storage.integerPartSize] != '.' && storage.integerPartSize < raw.size(); storage.integerPartSize++) {} 
    storage.integerPartSize -= !storage.isPositive; //If input string had a '-' sign we should start from 1 index
    //Without any loop we calculate size of fractional part by determining if the string has '.' and with using size of integer part
    storage.fractionalPartSize = raw.size() - storage.integerPartSize - (raw[storage.integerPartSize + !storage.isPositive] == '.') 
    - !storage.isPositive; 
    storage.integerPart = new int [storage.integerPartSize];
    storage.fractionalPart = new int [storage.fractionalPartSize];
    
    //If input string had a '-' sign we should start from 1 index
    for (int i = !storage.isPositive; i < storage.integerPartSize + !storage.isPositive; i++) { 
        storage.integerPart[i - !storage.isPositive] = raw[i] - '0'; // "- '0'" is a way to convert char digit to int type digit
    }

    for (int i = storage.integerPartSize + 1 + !storage.isPositive; i < raw.size(); i++) { //"i" is index of digit in the input string
        //"i - integerPartSize - 1" is index of digit in our dynamic array storage
        // "- '0'" is a way to convert char digit to int type digit
        storage.fractionalPart[i - storage.integerPartSize - 1 - !storage.isPositive] = raw[i] - '0'; 
    }

    return input; //Need for cascade using of the operator forces us to return input stream in the output of the operator function
}

ostream& operator<< (ostream& output, const BigNum& storage) {
    if (storage.isZero()) {
        cout << '0'; 
        return output;
    }

    if (!storage.isPositive) { //If the number was less than zero we want to print the minus sign before printing the number
        cout << '-';
    }
    
    bool flg = false; //It means we saw a non-zero digit in integerPart
    for (int i = 0; i < storage.integerPartSize; i++) {
        //If -1 was placed somewhere it means that the digit was considered for carry but currently is useless
        if (storage.integerPart[i] != -1 && (storage.integerPart[i] != 0 || flg || i == storage.integerPartSize - 1)) { 
            cout << storage.integerPart[i];
        }
        if (storage.integerPart[i] != 0) {
            flg = true;
        }
    }

    if (storage.fractionalPartSize != 0) { //Check if it is required to print '.' character in output
        cout << '.';
        for (int i = 0; i < min(storage.fractionalPartSize, MAXFRACSIZE); i++) { 
            cout << storage.fractionalPart[i];
        }
    }

    return output; //Need for cascade using of the operator forces us to return input stream in the output of the operator function
}

#endif