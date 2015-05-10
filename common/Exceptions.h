/**
@file
@brief Collection of Exceptions that can be thrown.
*/
#ifndef EXCEPTIONSH
#define EXCEPTIONSH

#include <iostream>
#include <string>

namespace tsquare
{

/**
@class
@brief Exception for attempting to access an out-of-bounds data element.
*/
class EOBException {

private:

std::string arrayname_;     /**< Name of the array accessed */
std::string classname_;     /**< Name of the class causing the exception */
std::string functionname_;  /**< Name of the function causing the exception */
int sizelimit_;        /**< Size of the array being accessed */
unsigned int indx_;    /**< Array element that was attempted */

public:

/**
@brief Default constructor.
*/
EOBException()
{
    classname_ = "";
    functionname_ = "";
    arrayname_ = "";
    sizelimit_ = 0;
    indx_ = 0;
}

/**
@brief Overloaded (usual) constructor.

Accepts parameters specifying what caused the exception.

@param &cname is the class name causing the exception
@param &fname is the function name causing the exception
@param &arname is the array name that went out of bounds
@param &sl is the size limit of the array
@param &id is the array element that was requested
*/
EOBException(const std::string &cname, const std::string &fname, const std::string &arname,
             const int sl, const unsigned int id)
{
    classname_ = cname;
    functionname_ = fname;
    arrayname_ = arname;
    sizelimit_ = sl;
    indx_ = id;
}

/**
@brief Gets the class name causing the exception.

@return the class name in which the out-of-bounds exception happened.
*/
std::string &getClassname() const { return (std::string &)classname_; }

/**
@brief Gets the function name causing the exception.

@return the function name in which the out-of-bounds exception happened.
*/
std::string &getFunctionname() const { return (std::string &)functionname_; }

/**
@brief Gets the array name that went out of bounds.

@return the array name that went out of bounds
*/
std::string &getArrayname() const { return (std::string &)arrayname_; }

/**
@brief Gets the size of the out-of-bounds array.

@return the size of the out-of-bounds array
*/
int getSizelimit() const { return sizelimit_; }

/**
@brief Gets the array element that was requested.

@return the requested array element
*/
unsigned int getIndx() const { return indx_; }

/**
@brief Prints a formatted copy of the exception details.
*/
void printException() {
    std::cout << std::endl << "EOB Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cerr << std::endl << "EOB Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    if (indx_ == 0) {
        std::cout << "        Array: " << arrayname_ << std::endl;
        std::cerr << "        Array: " << arrayname_ << std::endl;
    } else {
        std::cout << "        Array: " << arrayname_ << " contains " << sizelimit_;
        std::cout << " elements, but tried to access element " << indx_ << std::endl;
        std::cerr << "        Array: " << arrayname_ << " contains " << sizelimit_;
        std::cerr << " elements, but tried to access element " << indx_ << std::endl;
    }
    return;
}
};


/**
@class
@brief Exception for file I/O errors.
*/
class FileException {

private:

std::string filename_;      /**< Name of the file causing the exception */
std::string classname_;     /**< Name of the class causing the exception */
std::string functionname_;  /**< Name of the function causing the exception */
std::string extype_;        /**< Type of exception */

public:

/**
@brief Default constructor.
*/
FileException() {
    classname_ = "";
    functionname_ = "";
    filename_ = "";
    extype_ = "";
}

/**
@brief Overloaded (usual) constructor.

Accepts parameters specifying what caused the exception.

@param &cname is the class name causing the exception
@param &fname is the function name causing the exception
@param &filename is the file that caused the exception
@param &extype is the exception type
*/
FileException(const std::string &cname, const std::string &fname,
              const std::string &filename, const std::string &extype) {
    classname_ = cname;
    functionname_ = fname;
    filename_ = filename;
    extype_ = extype;
}

/**
@brief Gets the class name causing the exception.

@return the class name in which the out-of-bounds exception happened.
*/
std::string &getClassname() const { return (std::string &)classname_; }

/**
@brief Gets the function name causing the exception.

@return the function name in which the out-of-bounds exception happened.
*/
std::string &getFunctionname() const { return (std::string &)functionname_; }

/**
@brief Gets the file causing the exception.

@return the file in which the out-of-bounds exception happened.
*/
std::string &getFilename() const { return (std::string &)filename_; }

/**
@brief Gets the description of the exception.

@return a string description of the exception
*/
std::string &getExtype() const { return (std::string &)extype_; }

/**
@brief Prints a formatted copy of the exception details.
*/
void printException() {
    std::cout << std::endl << "File Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cout << "        File: " << filename_ << ", Problem:" << extype_ << std::endl;
    std::cerr << std::endl << "GEM Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cerr << "        File: " << filename_ << ", Problem: " << extype_ << std::endl;
    return;
}
};

/**
@class
@brief Exception for floating point errors.
*/
class FloatException {

private:

std::string description_;      /**< Description of the floating-point exception */
std::string classname_;        /**< Name of the class causing the exception */
std::string functionname_;     /**< Name of the function causing the exception */

public:

/**
@brief Default constructor.
*/
FloatException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
}

/**
@brief Overloaded (usual) constructor.

Accepts parameters specifying what caused the exception.

@param &cname is the class name causing the exception
@param &fname is the function name causing the exception
@param &strd is the string description of the exception
*/
FloatException(const std::string &cname, const std::string &fname, const std::string &strd) {
    classname_ = cname;
    functionname_ = fname;
    description_ = strd;
}

/**
@brief Gets the class name causing the exception.

@return the class name in which the out-of-bounds exception happened.
*/
std::string &getClassname() const { return (std::string &)classname_; }

/**
@brief Gets the function name causing the exception.

@return the function name in which the out-of-bounds exception happened.
*/
std::string &getFunctionname() const { return (std::string &)functionname_; }

/**
@brief Gets the description of the exception.

@return a string description of the exception
*/
std::string &getDescription() const { return (std::string &)description_; }

/**
@brief Prints a formatted copy of the exception details.
*/
void printException() {
    std::cout << std::endl << "Floating Point Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cout << "        Description: " << description_ << std::endl;
    std::cerr << std::endl << "Floating Point Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cerr << "        Description: " << description_ << std::endl;
    return;
}
};

/**
@class
@brief Exception for bad data of one kind or another.
*/
class DataException {

private:

std::string description_;      /**< Description of the data exception */
std::string classname_;        /**< Name of the class causing the exception */
std::string functionname_;     /**< Name of the function causing the exception */

public:

/**
@brief Default constructor.
*/
DataException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
}

/**
@brief Overloaded (usual) constructor.

Accepts parameters specifying what caused the exception.

@param &cname is the class name causing the exception
@param &fname is the function name causing the exception
@param &strd is the string description of the exception
*/
DataException(const std::string &cname, const std::string &fname, const std::string &strd) {
    classname_ = cname;
    functionname_ = fname;
    description_ = strd;
}

    // Getters

/**
@brief Gets the class name causing the exception.

@return the class name in which the out-of-bounds exception happened.
*/
std::string &getClassname() const { return (std::string &)classname_; }

/**
@brief Gets the function name causing the exception.

@return the function name in which the out-of-bounds exception happened.
*/
std::string &getFunctionname() const { return (std::string &)functionname_; }

/**
@brief Gets the description of the exception.

@return a string description of the exception
*/
std::string &getDescription() const { return (std::string &)description_; }

/**
@brief Prints a formatted copy of the exception details.
*/
void printException() {
    std::cout << std::endl << "Data Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cout << "        Problem:" << description_ << std::endl;
    std::cerr << std::endl << "Data Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_ << "::" << functionname_ << std::endl;
    std::cerr << "        Problem: " << description_ << std::endl;
    return;
}
};

} // End of namespace tsquare
#endif
