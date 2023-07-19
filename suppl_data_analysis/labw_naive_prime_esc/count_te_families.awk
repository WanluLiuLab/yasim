BEGIN {
    OFS="\t"
}
{
    count[$1]++
}
END {
    for (word in count){
        print word, count[word]
    }
}
