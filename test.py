import os
import platform

os_type = platform.system()
spacia_path = os.path.abspath(__file__)
spacia_fn = spacia_path.replace('test.py', 'spacia.py')
input_path = spacia_path.replace('test.py', 'test/input')
test_path = spacia_path.replace('test.py', 'test')
counts_fn = os.path.join(input_path,'counts.txt')
meta_fn = os.path.join(input_path,'spacia_metadata.txt')

print('Testing Spacia with a single gene as response feature and simple aggregation')
params = '-rc A -sc B -rf gene1 -sf gene2,gene3 -d 5 -m 2000,1000,10,1 -nc 20'
output_path = os.path.join(test_path,'single_gene_simple_agg')
cmd = 'python {} {} {} {} -o {}'.format(
    spacia_fn, counts_fn, meta_fn, params, output_path)
codes = os.system(cmd)
if codes == 0:
    print('Test Succeeded.')
else:
    print('Test failed, please check log at {}'.format(output_path))

print('Testing Spacia with multiple genes as response feature and no agg mode')
output_path = os.path.join(test_path,'multi_gene_no_agg')
params = '-rc A -sc B -rf "gene1|gene2|gene3" -sf gene2,gene3 -d 5 -m 2000,1000,10,1 -rec auto'
cmd = 'python {} {} {} {} -o {}'.format(
    spacia_fn, counts_fn, meta_fn, params, output_path)
codes = os.system(cmd)
if codes == 0:
    print('Test Succeeded.')
else:
    print('Test failed, please check log at {}'.format(output_path))

print('Testing Spacia with multiple genes as response feature and pca agg mode')
output_path = os.path.join(test_path,'multi_gene_pca_agg')
params = '-rc A -sc B -rf "gene1|gene2|gene3" -sf pca -d 5 -m 2000,1000,10,1'
cmd = 'python {} {} {} {} -o {}'.format(
    spacia_fn, counts_fn, meta_fn, params, output_path)
codes = os.system(cmd)
if codes == 0:
    print('Test Succeeded.')
else:
    print('Test failed, please check log at {}'.format(output_path))
        
