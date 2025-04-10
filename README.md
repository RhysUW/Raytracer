breif explanation: 
First program parses the input via parseFile() to load all the scene parameters . Each shperes transformation matrix is build in buildSphereMatrix() by combining a scale and translation, and its inverse is computed for accurate ray-shpere intersection and normal calculation. in traceRay() we check each sphere for intersection (via intersectSphere()), then compute local illumination using ambient, diffuse, specular and clamp the results to [0...1] 

iv only done two of the tests, i really dont even want to see if the other ones fail i might drive into oncoming traffic.

theres still a bunch of printf statments i used for debugging, to lazy to remove them sorry.