from netCDF4 import Dataset

new_data = Dataset('/net/home/h01/tbendall/results/vector_transport_paper/conv_1_quads_cylinder/global_output.nc', 'w')
new_data.createDimension('run_id', None)
new_data.createDimension('time', 1)
new_data.createVariable('run_id', int, ('run_id',))
new_data.createVariable('time', float, ('time',))
new_data.createGroup('F_0')
new_data.createGroup('F_1')
new_data.createGroup('F_2')
new_data['F_0'].createGroup('errors')
new_data['F_1'].createGroup('errors')
new_data['F_2'].createGroup('errors')
new_data['F_0']['errors'].createVariable('L2_error', float, ('run_id', 'time'))
new_data['F_1']['errors'].createVariable('L2_error', float, ('run_id', 'time'))
new_data['F_2']['errors'].createVariable('L2_error', float, ('run_id', 'time'))

new_data.createVariable('dx', float, ('run_id',))
new_data.createVariable('dt', float, ('run_id',))

new_data['time'][0] = 100.0

new_data.close()

old_files = ['a', 'b', 'c']

new_data = Dataset('/net/home/h01/tbendall/results/vector_transport_paper/conv_1_quads_cylinder/global_output.nc', 'a')

for old_file in old_files:
    old_data = Dataset(f'/net/home/h01/tbendall/results/vector_transport_paper/conv_1{old_file}_cyl/global_output.nc', 'r')
    
    for i, old_run_id in enumerate(old_data['run_id'][:]):
        new_run_id = len(new_data['run_id'][:])
        new_data['run_id'][new_run_id] = new_run_id
        new_data['dx'][new_run_id] = old_data['dx'][i]
        new_data['F_0']['errors']['L2_error'][new_run_id, -1] = old_data['F_0']['errors']['L2_error'][i,-1]
        new_data['F_1']['errors']['L2_error'][new_run_id, -1] = old_data['F_1']['errors']['L2_error'][i,-1]
        new_data['F_2']['errors']['L2_error'][new_run_id, -1] = old_data['F_2']['errors']['L2_error'][i,-1]

    old_data.close()

new_data.close()