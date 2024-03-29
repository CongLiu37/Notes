#!/usr/bin/env python3

import sys, os

ds = os.listdir(sys.argv[1])

for d in ds:
    fPath = sys.argv[1] + '/' + d + '/evm.out'
    size = os.path.getsize(fPath)
    if size > 0:
        blocks = open(fPath).read().strip().split('#')[1:]
        for block in blocks:
            coords = []
            evidence = []
            for line in block.strip().split('\n')[1:]:
                if line.strip() != '' and line[0] != '!':
                    meta = line.strip().split('\t')
                    coords.append(int(meta[0]))
                    coords.append(int(meta[1]))
                    coords.sort()
                    evidence.extend([tuple(x[1:-1].split(';')) for x in meta[-1].split(',')])

            evidence = set(evidence)
            sources = set([x[1] for x in evidence])

            print(d + '\t' + 'EVM' + '\t'+ 'gene' + '\t' + str(coords[0]) + '\t' + str(coords[-1]) + '\t' + ','.join([x[0] for x in evidence]) + '\t' + ','.join(sources))
            
# #!/usr/bin/env python3
# 
# import sys, os
# 
# ds = os.listdir(sys.argv[1])
# 
# for d in ds:
#     fPath = sys.argv[1] + '/' + d + '/evm.out'
#     if os.path.exists(fPath):
#         size = os.path.getsize(fPath)
#         if size > 0:
#             os.system("sed 's/##/\*\*/' "+fPath+" > "+fPath+".edit")
#             fPath=fPath+".edit"
#             blocks = open(fPath).read().strip().split('#')[1:]
#             for block in blocks:
#                 coords = []
#                 evidence = []
#                 for line in block.strip().split('\n')[1:]:
#                     if line.strip() != '' and line[0] != '!':
#                         meta = line.strip().split('\t')
#                         coords.append(int(meta[0]))
#                         coords.append(int(meta[1]))
#                         coords.sort()
#                         ev=[]
#                         for x in meta[-1].split(','):
#                             if x[-1]=='}' and x[0]=='{':
#                                 ev.append(x)
#                 evidence.extend([tuple(x[1:-1].split(';')) for x in ev])
#                       #evidence.extend([tuple(x[1:-1].split(';')) for x in meta[-1].split(',')])
# 
#                 evidence = set(evidence)
#                 sources = set([x[1] for x in evidence])
#                 print(d + '\t' + 'EVM' + '\t'+ 'gene' + '\t' + str(coords[0]) + '\t' + str(coords[-1]) + '\t' + ','.join([x[0] for x in evidence]) + '\t' + ','.join(sources))
#                 #if coords!=[] and evidence!=[]:
#                 #    print(d + '\t' + 'EVM' + '\t'+ 'gene' + '\t' + str(coords[0]) + '\t' + str(coords[-1]) + '\t' + ','.join([x[0] for x in evidence]) + '\t' + ','.join(sources))
#             
#               
# # fPath="temp_dir//contig_10132/evm.out"
# # blocks = open(fPath).read().strip().split('#')[1:]
# # for block in blocks:
# #     coords = []
# #     evidence = []
# #     for line in block.strip().split('\n')[1:]:
# #         if line.strip() != '' and line[0] != '!':
# #             meta = line.strip().split('\t')
# #             coords.append(int(meta[0]))
# #             coords.append(int(meta[1]))
# #             coords.sort()
# #             ev=[]
# #             for x in meta[-1].split(','):
# #                 if x[-1]=='}' and x[0]=='{':
# #                     ev.append(x)
# #     evidence.extend([tuple(x[1:-1].split(';')) for x in ev])
# #           #evidence.extend([tuple(x[1:-1].split(';')) for x in meta[-1].split(',')])
# #     evidence = set(evidence)
# #     sources = set([x[1] for x in evidence])
# #     print(fPath)
# #     print(coords)
