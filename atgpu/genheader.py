#Generate include file for GPU generated code
print("Generating atgpu/element.gpuh")

#Read the include file
fr = open("atgpu/Element.h","r")
lines = fr.readlines()
fr.close()

# Add type definition
idx = [idx for idx, element in enumerate(lines) if element.startswith('// DEFINE_TYPE')]
if len(idx)==0:
    raise RuntimeError("Invalid atgpu/Element.h file \'// DEFINE_TYPE\' not found")
I = idx[0]
lines.pop(I)
lines.insert(I,"typedef signed char        int8_t;\n")
lines.insert(I,"typedef signed short       int16_t;\n")
lines.insert(I,"typedef signed int         int32_t;\n")
lines.insert(I,"typedef signed long long   int64_t;\n")
lines.insert(I,"typedef unsigned char      uint8_t;\n")
lines.insert(I,"typedef unsigned short     uint16_t;\n")
lines.insert(I,"typedef unsigned int       uint32_t;\n")
lines.insert(I,"typedef unsigned long long uint64_t;\n")

# Remove comments
lines = [l for l in lines if not l.startswith('//')]

#insert raw string marker
lines.insert(0,"R\"(\n")
lines.insert(len(lines),"\n)\"\n")

# Write the file
fw = open("atgpu/element.gpuh","w")
fw.writelines(lines)
fw.close()
