#Find all files with extention ".java"

FILES=$(find src -type f -name '*.java')

#Iterate over all the files

for i in $FILES
    do 
      echo "cleaning file " $i
      #sed -i '' '/^@TestClass/d' $i
      #sed -i '' '/^@TestMethod/d' $i
      sed -i '' '/^.*TestClass.*$/d' $i
      sed -i '' '/^.*TestMethod.*$/d' $i
    done

#sed -i '' '/^@TestClass/d' org/openscience/smsd/*.java
    
