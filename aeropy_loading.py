import aeropy.xfoil_module as xf

def get_pressure_distribution(foil):
    '''
    oh my this is so stupid, that you have to do this is maddening and took me
    so long to figure out i just when and wrote this function
    and completly went around how your meant to do it
    if i could re do my dissertation i would just fix this disaster by rewriting xfoil in python
    
    dont touch anything or it will 100% break
    '''
    foil_name = "do_not_delete_this_or_everything_breaks"
    for line in foil:
        output_file.write("     "+str(format(line[0],'.6f'))+"    "+str(format(line[1],'.6f'))+"\n")#
    output_file.close()
    Data = xf.find_pressure_coefficients(foil_name, 0., iteration=100, NACA=False)
    return Data
