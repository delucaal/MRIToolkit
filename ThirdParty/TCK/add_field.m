% Taken from MRtrix
function image = add_field (image, key, value)
  if isfield (image, key)
    previous = getfield (image, key);
    if iscell (previous)
      image = setfield (image, key, [ previous {value} ]);
    else
      image = setfield (image, key, { previous, value });
    end
  else
    image = setfield (image, key, value);
  end
end