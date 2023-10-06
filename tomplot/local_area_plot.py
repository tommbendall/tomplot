import pandas as pd

#TODO
# Perhaps instead of an area function we could just call this plotting area and 
# if no limits are passed it applied the gusto / lfric grid
def area_restriction(field_data, coords_X, coords_Y,  Y_lim=None, X_lim=None):
    """
    a function that takes a data dictionary and will restrict it for a given value
    For the coordinate fields you must pass in a 1d object so either coords[:,level]
    for slicing on z or coords.flatten if slicing along lat / lon

    args:
    field_data: The Plotting field you wish to return 
    coords_X: The coordinates you wish to plot along the X axis.
    coords_Y: The field you wish to plot along the Y axis.
    X_lim: tuple containing the (lower, upper) limits for the X axis
    Y_lim: tuple containing the (lower, upper) limits for the Y axis
    returns
    restricted field data, X_coords, Y_coords
    """
    data_dict = {'field': field_data, 'X': coords_X, 'Y': coords_Y}
    df = pd.DataFrame(data_dict) 
    # ive used this logic to allow a user to only enter a limit for one axis
    # and not have to use a new function, im sure there is a more elegant way
    if not X_lim is None:
        X_min, X_max = X_lim
        df = df[(df["X"] >= X_min) & (df["X"] <= X_max)]

    if not Y_lim is None:
        Y_min, Y_max = Y_lim
        df = df[(df["Y"] >= Y_min) & (df["Y"] <= Y_max)]
        
    new_field_data = df['field']
    new_coords_X = df['X']
    new_coords_Y = df['Y']

    return (new_field_data, new_coords_X, new_coords_Y)

    
    