import numpy as np

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()


d=open('GSE81682_HTSeq_counts_header.txt', 'r')

d0 = d.readline().split()

label = []
ct_id_dict = {}
ct_id_mat = []
i = 0
for ct_id in d0[1:]:
	ct = ct_id.split('_')[0]
	if ct in ct_id_dict:
		label.append(ct_id_dict[ct])
	else:
		i = i+1
		ct_id_dict[ct] = i
		ct_id_mat.append([ct, i])
		label.append(i)

label = np.array(label)
label = label.reshape((label.shape[0],1))
ct_id_mat = np.array(ct_id_mat)


write2d_array(label, 'GSE81682_label.txt')
write2d_array(ct_id_mat, 'GSE81682_ct2label.txt')



