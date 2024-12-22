template <typename T>
T abs(T const a){
    return a > 0 ? a : -a;
}

template <typename T>
const T abs(T const  a)
{
    return a > 0 ? a : -a;
}