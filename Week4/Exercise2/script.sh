for step in 8 16 32 64
do
cat > input << EOF
500
$step
0
0
0
0
0
0
0
0
0
0
0
0
EOF
./a.out < input >> dati
done

