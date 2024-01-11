"all vars to type T"
convert_many(T::Type, vars...) = Tuple(convert(T, x) for x in vars)

"promoted type of all vars"
promote_type_of_many(vars...) = promote_type(Tuple(typeof(x) for x in vars)...)

