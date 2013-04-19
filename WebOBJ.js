/******

Infect! loadObj.js

Description:
Rewrite of loadObj.h/loadObj.cpp for class.
Allows user to load in a .obj file if one exists.

@author: Sean Lander

******/

var objm = objm || { };

objm.C_PI = 3.14159265;

objm.OBJM_NONE =         0;   // render with vertices
objm.OBJM_FLAT_SHADE =   1;   // render with face normals
objm.OBJM_SMOOTH_SHADE = 2;   // render with vertex normals
objm.OBJM_TEXTURE =      4;   // render with texture coords
objm.OBJM_COLOR =        8;   // render with colors
objm.OBJM_MATERIAL =     16;  // render with materials

/*
Easy access function to get a model m's triangle x's attributes,
such as vindices, nindices, tindices
*/
objm.T = function(m,x,a) { return m.triangles[x][a].elements; };
    
/*
OBJMmaterial: Structure that defines a material in a model.
Read in from the appropriate line in the OBJ file, this structure
contains info that defines how the object appears when lighted.
*/
objm.OBJMmaterial = function()
{
    this.name =      "";             // Name of the material
    this.ambient =   Vector.Zero(4); // Ambient component
    this.diffuse =   Vector.Zero(4); // Diffuse component
    this.specular =  Vector.Zero(4); // Specular component
    this.emissive =  Vector.Zero(4); // Emissive component
    this.shininess = 0;              // Specular exponent
};

/*
OBJMtriangle: Structure that defines a triangle in a model. All
the information relevant to a triangle, which is the building block
of any OBJ mesh. Array of these is used in the OBJMmodel structure.
*/
objm.OBJMtriangle = function()
{
    this.vindices = Vector.Zero(3); // Array of triangle vertex indices
    this.nindices = Vector.Zero(3); // Array of triangle normal indices
    this.tindices = Vector.Zero(3); // Array of triangle texcoord indices
    this.findex =   0;              // Triangle facet normal indices
};

/*
OBJMgroup: Strcuture that defines a group in a model. Large OBJ files
representing complex meshes have groups in theme, which is a means of
grouping different regions of the mesh. For example, an OBJ mesh of a
human might have the triangles for the head, torso, arms and legs in
different groups
*/
objm.OBJMgroup = function()
{
    this.name =      "";   // The group name
    this.triangles = [];   // GLuint*    // The array of triangle indices
    this.material =  0;    // GLuint     // The index of the material for the group
    this.next =      null; // OBJMgroup* // The next group in the model
    
    this.numtriangles = function() { return this.triangles.length; } // GLuint // The number of triangles in the group
};

/*
OBJMmodel: Structure that defines a model. It has arrays for all the info
needed when loading an OBJ file:

Vertices
Face normals
Vertex normals
Texture coordinates
Triangles
Materials

It also has a linked list for the groups in the OBJ file.
*/
objm.OBJMmodel = function()
{
    this.pathname =    ""; // char* // Path to the model
    this.mtllibname =  ""; // char* // Name of the material library
    
    this.vertices =    []; // GLfloat*      // Array of vertices
    this.facetnorms =  []; // GLfloat*      // Array of facetnorms
    this.vertexnorms = []; // GLfloat*      // Array of vertexnorms
    this.texcoords =   []; // GLfloat*      // Array of texcoords
    this.triangles =   []; // OBJMtriangle* // Array of traingles
    this.materials =   []; // OBJMmaterial* // Array of materials
    this.groups =      []; // OBJMgroup*    // Linked list of groups

    this.numvertices =    function() { return this.vertices.length; };    // GLuint // Number of vertices
    this.numfacestnorms = function() { return this.facetnorms.length; };  // GLuint // Number of facetnorms
    this.numvertexnorms = function() { return this.vertexnorms.length; }; // GLuint // Number of vertexnorms
    this.numtexcoords =   function() { return this.texcoords.length; };   // GLuint // Number of texcoords in the model
    this.numtriangles =   function() { return this.triangles.length; };   // GLuint // Number of triangles in the model
    this.nummaterials =   function() { return this.materials.length; };   // GLuint // Number of materials in the model
    this.numgroups =      function() { return this.groups.length; };      // GLuint // Number of groups in the model

    this.position = Vector.Zero(3); // GLfloat // Position of the model
};

/*
OBJMnode: General purpose node
*/
objm.OBJMnode = function()
{
    this.index =    0;     // Index of node
    this.averaged = false; // ???
    this.next =     null;  // Next node
}
objm.loadObj = function() {
    /* Public */
    this.eye =          Vector.Zero(3); // GLfloat[3]
    this.at =           Vector.Zero(3); // GLfloat[3]
    this.up =           Vector.Zero(3); // GLfloat[3]
    this.translation =  Vector.Zero(3); // GLfloat[3]
    this.rotation =     Vector.Zero(4); // GLfloat[4]

    this.projection =   Vector.Zero(16); // GLdouble[16]
    this.modelview =    Vector.Zero(16); // GLdouble[16]
    this.inverse =      Vector.Zero(16); // GLdouble[16]

    this.swapped =      false; // GLboolean
    this.world_draw =   false; // GLboolean
    this.pmodel =       null;  // OBJMmodel
    this.selection =    0;     // GLint

    /* Private */
    var xRot = 0; // int
    var yRoy = 0; // int
    var zRot = 0; // int
    
    var drawType =          2;    // int
    var drawPointsAlso =    0;    // int
    var nearPlane =         0.1;  // float (incorrect spelling in source)
    var farPlane =          20;   // float

    var curPickDepth =      0;              // GLfloat
    var minVertIndex =      0;              // GLint
    var minVertCoords =     Vector.Zero(3); // GLfloat[3]
    var distClosestVertex = 0;              // GLfloat
    
    var viewport =      Vector.Zero(4);  // GLint[4]
    var mvmatrix =      Vector.Zero(16); // GLdouble[16]
    var projmatrix =    Vector.Zero(16); // GLdouble[16]
    var edgeList =      null;            // vector< vector<int> >
    // TODO: Need to make vector class
}

/*
Scales the model to fit within a cube of dimensions 2 x 2 x 2. This way, the model is
centered at (0, 0, 0) and the 'imaginary bounding cube' has corners at (1, 1, 1) and
at (-1, -1, -1). This is necessary if a single scene and camera model are to be used
for any loaded object, since objects have different sizes.
*/
objm.loadObj.prototype.objmUnitize = function(model) { // (OBJMmodel*)
    var i; // GLuint
    var maxx, minx, maxy, miny, maxz, minz; // GLfloat
    var cx, cy, cz, w, h, d; // GLfloat
    var scale; // GLfloat

    if(assert(model, 'Model exists')) return;
    if(assert(model.vertices, 'Model has vertices')) return;

    // Get the maximums/minimums
    maxx = minx = model.vertices[3 + 0];
    maxy = miny = model.vertices[3 + 1];
    maxz = minz = model.vertices[3 + 2];

    for(i = 1 ; i <= model.numvertices() ; i++) {
        maxx = (maxx < model.vertices[3 * i + 0] ? model.vertices[3 * i + 0] : maxx);
        minx = (minx > model.vertices[3 * i + 0] ? model.vertices[3 * i + 0] : maxx);

        maxy = (maxy < model.vertices[3 * i + 1] ? model.vertices[3 * i + 1] : maxy);
        miny = (miny > model.vertices[3 * i + 1] ? model.vertices[3 * i + 1] : miny);

        maxz = (maxz < model.vertices[3 * i + 2] ? model.vertices[3 * i + 2] : maxz);
        minz = (minz > model.vertices[3 * i + 2] ? model.vertices[3 * i + 2] : minz);
    }

    // Calculate the model width, height, and depth
    w = this.objmAbs(maxx) + this.objmAbs(minx);
    h = this.objmAbs(maxy) + this.objmAbs(miny);
    d = this.objmAbs(maxz) + this.objmAbs(minz);

    // Calculate the center of the model
    cx = (maxx + minx) / 2;
    cy = (maxy + miny) / 2;
    cz = (maxz + minz) / 2;

    // Calculate the unitizing scale factor
    scale = 2.0 / this.objmMax(this.objmMax(w, h), d);

    // Translate to the center then scale
    for(i = 1 ; i <= model.numvertices() ; i++) {
        model.vertices[3 * i + 0] -= cx;
        model.vertices[3 * i + 1] -= cy;
        model.vertices[3 * i + 2] -= cz;
        model.vertices[3 * i + 0] *= scale;
        model.vertices[3 * i + 1] *= scale;
        model.vertices[3 * i + 2] *= scale;
    }

    return model;
}; // GLfloat

/*
Gets the dimensions of the object - l, b, h. The final values
in the dimensions array represents the bounding box for the object.
*/
objm.loadObj.prototype.objmDimensions = function(model, dimensions) { // (OBJMmodel*, GLfloat*)
    var i; // GLuint
    var maxx, minx, maxy, miny, maxz, minz; // GLfloat

    if(assert(model, 'Model exists')) return;
    if(assert(model.vertices, 'Model has vertices')) return;
    if(assert(dimensions, 'Dimensions exist')) return;

    // Get the maximums/minimums
    maxx = minx = model.vertices[3 + 0];
    maxy = miny = model.vertices[3 + 1];
    maxz = minz = model.vertices[3 + 2];

    for(i = 1 ; i <= model.numvertices() ; i++) {
        maxx = (maxx < model.vertices[3 * i + 0] ? model.vertices[3 * i + 0] : maxx);
        minx = (minx > model.vertices[3 * i + 0] ? model.vertices[3 * i + 0] : maxx);

        maxy = (maxy < model.vertices[3 * i + 1] ? model.vertices[3 * i + 1] : maxy);
        miny = (miny > model.vertices[3 * i + 1] ? model.vertices[3 * i + 1] : miny);

        maxz = (maxz < model.vertices[3 * i + 2] ? model.vertices[3 * i + 2] : maxz);
        minz = (minz > model.vertices[3 * i + 2] ? model.vertices[3 * i + 2] : minz);
    }

    // Calculate the model width, height, and depth
    dimensions[0] = this.objmAbs(maxx) + this.objmAbs(minx);
    dimensions[1] = this.objmAbs(maxy) + this.objmAbs(miny);
    dimensions[2] = this.objmAbs(maxz) + this.objmAbs(minz);

}; // GLvoid

/*
Scales the model by the passed in scale factor
*/
objm.loadObj.prototype.objmScale = function(model, scale) { // (OBJMmodel*, GLfloat)
    var i; // GLuint

    for(i = 1 ; i < model.numvertices() ; i++) {
        model.vertices[3 * i + 0] *= scale;
        model.vertices[3 * i + 1] *= scale;
        model.vertices[3 * i + 2] *= scale;
    }
}; // GLvoid

/*
Reverses the winding of the triangle faces of the model. This
is done by changing the vertex order - 1,2,3 becomes 3,2,1; this
effectively flips the faces when CCW orientation is used. The normals
are pointed in the opposite direction. Thus, the net effect is that
the model is turned 'inside out'.
*/
objm.loadObj.prototype.objmReverseWinding = function(model) { // (OBJMmodel*)
    var i, swap; // GLuint

    if(assert(model, 'Model exists')) return;

    // Reverse the winding of the faces
    for(i = 0 ; i < model.numtriangles() ; i++) {
        swap = objm.T(model,i,'vindices')[2];
        objm.T(model,i,'vindices')[0] = objm.T(model,i,'vindices')[2];
        objm.T(model,i,'vindices')[2] = swap;

        if(model.numvertexnorms()) {
            swap = objm.T(model,i,'nindices')[0];
            objm.T(model,i,'nindices')[0] = objm.T(model,i,'nindices')[2];
            objm.T(model,i,'nindices')[2] = swap;
        }

        if(model.numtexcoords()) {
            swap = objm.T(model,i,'tindices')[0];
            objm.T(model,i,'tindices')[0] = objm.T(model,i,'tindices')[2];
            objm.T(model,i,'tindices')[2] = swap;
        }
    }

    // Reverse the facet normals
    for(i = 1 ; i <= model.numfacestnorms() ; i++) {
        model.facetnorms[3 * i + 0] *= -1;
        model.facetnorms[3 * i + 1] *= -1;
        model.facetnorms[3 * i + 2] *= -1;
    }

    // Reverse the vertex normals
    for(i = 1 ; i <= model.numvertexnorms() ; i++) {
        model.vertexnorms[3 * i + 0] *= -1;
        model.vertexnorms[3 * i + 1] *= -1;
        model.vertexnorms[3 * i + 2] *= -1;
    }

}; // GLvoid

/*
Get the face normals from the OBJ file and normalize them.
Normalization is necessary when lighting needs to be done on
the model;
*/
objm.loadObj.prototype.objmFacetNormals = function(model) { // (OBJMmodel*)
    var i; // GLuint
    u = Vector.Zero(3); // GLfloat[3]
    v = Vector.Zero(3); // Glfloat[3]

    if(assert(model, 'Model exists')) return;
    if(assert(model.vertices, 'Model vertices exists')) return;

    // Free any old facetnormals
    // Allocatte memory for the new facet normals
    // This is Javascript! Screw memory management!

    for(i = 0 ; i < model.numtriangles() ; i++) {
        model.triangles[i].findex = i+1;

        u.elements[0] = model.vertices[3 * objm.T(model,i,'vindices')[1] + 0] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 0];
        u.elements[1] = model.vertices[3 * objm.T(model,i,'vindices')[1] + 1] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 1];
        u.elements[2] = model.vertices[3 * objm.T(model,i,'vindices')[1] + 2] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 2];

        v.elements[0] = model.vertices[3 * objm.T(model,i,'vindices')[2] + 0] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 0];
        v.elements[1] = model.vertices[3 * objm.T(model,i,'vindices')[2] + 1] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 1];
        v.elements[2] = model.vertices[3 * objm.T(model,i,'vindices')[2] + 2] - model.vertices[3 * objm.T(model,i,'vindices')[0] + 2];

        this.objmCross(u, v, model.facetnorms[3 * (i+1)]);
        this.objmNormalize(model.facetnorms[3 * (i+1)]);
    }
}; // GLvoid

// Computes the normals at every vertex.
objm.loadObj.prototype.objmVertexNormals = function(model, angle) { // (OBJMmodel*, GLfloat)
    var node; // OBJMnode*
    var tail; // OBJMnode*
    var members = []; // OBJMnode**
    var normals = []; // GLfloat*
    var numnormals = function() { return normals.length; }; // GLuint
    var average = [0,0,0]; // GLfloat[3]
    var dot, cos_angle; // GLfloat
    var i, avg; // GLuint

    if(assert(model, 'Model exists')) return;
    if(assert(model.facetnorms, 'Model has face normals')) return;

    // Calculate the cosine of the angle (in degrees)
    cos_angle = Math.cos(angle * C_PI / 180);

    // Free any previous normals
    // Allocate space for the new normals
    // Allocate a structure that will hold a linked list of triangle indicies for each vertex
    // This is Javascript! Screw memory management!

    // Linked list of triangle indices for each vertex
    members.push(null);
    for(i = 1 ; i < model.numvertices() ; i++) {
        members.push(null);
    }

    for(i = 0 ; i < model.numtriangles() ; i++) {
        node = new objm.OBJMnode();
        node.index = i;
        node.next = members[objm.T(model,i,'vindices')[0]];
        members[objm.T(model,i,'vindices')[0]] = node;

        node = new objm.OBJMnode();
        node.index = i;
        node.next = members[objm.T(model,i,'vindices')[1]];
        members[objm.T(model,i,'vindices')[1]] = node;

        node = new objm.OBJMnode();
        node.index = i;
        node.next = members[objm.T(model,i,'vindices')[2]];
        members[objm.T(model,i,'vindices')[2]] = node;
    }

    // Calculate the average normal for each vertex
    numnormals = 1;
    for(i = 1 ; i <= model.numvertices() ; i++) {
        // Calculate an average normal for this vertex by averaging the
        // facet normal of every triangle this vertex is in
        node = members[i];

        if(!node) console.error("objmVertexNormals(): vertex w/o a triangle");

        average[0] = 0; average[1] = 0; average[2] = 0;
        avg = 0;

        while(node) {
            // only average if the dot product of the angle between the two
            // facet normals is greater than the cosine of the threshold
            // angle -- or, said another way, the angle between the two
            // facet normals is less than (or equal to) the threshold angle

            dot = this.objmDot(model.facetnorms[3 * model.triangles[node.index].findex, model.facetnorms[3 * model.triangles[members[i].index].findex]);
            if(dot > cos_angle) {
                node.averaged = true;
                average[0] += model.facetnorms[3 * model.triangles[node.index].findex + 0];
                average[1] += model.facetnorms[3 * model.triangles[node.index].findex + 1];
                average[2] += model.facetnorms[3 * model.triangles[node.index].findex + 2];
            
                avg = 1; // At least one normal averaged
            } else {
                node.averaged = false;
            }
            node = node.next;
        }

        if(avg) {
            // normalize the averaged normal
            this.objmNormalize(average);

            // add the normal to the vertex normals list
            model.vertexnorms[3 * numnormals + 0] = average[0];
            model.vertexnorms[3 * numnormals + 1] = average[1];
            model.vertexnorms[3 * numnormals + 2] = average[2];
            avg = numnormals;
            numnormals++;
        }

        // Set the normal of this vertex in each triangle it is in
        node = members[i];

        while(node) {
            if(node.averaged) {
                // If this node was averaged, use the average normal
                if(objm.T(model,node.index,'vindices')[0] == i)
                    objm.T(model,node.index,'vindices')[0] = avg;
                else if(objm.T(model,node.index,'vindices')[1] == i)
                    objm.T(model,node.index,'vindices')[0] = avg;
                else if(objm.T(model,node.index,'vindices')[2] == i)
                    objm.T(model,node.index,'vindices')[2] = avg;
            } else {
                model.vertexnorms[3 * numnormals + 0] = model.facetnorms[3 * model.triangles[node.index].findex + 0];
                model.vertexnorms[3 * numnormals + 1] = model.facetnorms[3 * model.triangles[node.index].findex + 1];
                model.vertexnorms[3 * numnormals + 2] = model.facetnorms[3 * model.triangles[node.index].findex + 2];

                if(objm.T(model,node.index,'vindices')[0] == i)
                    objm.T(model,node.index,'vindices')[0] = numnormals;
                else if(objm.T(model,node.index,'vindices')[1] == i)
                    objm.T(model,node.index,'vindices')[0] = numnormals;
                else if(objm.T(model,node.index,'vindices')[2] == i)
                    objm.T(model,node.index,'vindices')[2] = numnormals;
                numnormals++;
            }
            node = node.next;
        }
    }

    model.numvertexnorms = numnormals - 1;

    // Free the member information
    // Psha

    // Pack the normals array (we previously allocated the maximum
    // number of normals that could possibly be created (numtriangles *
    // 3), so get rid of some of them (usually a lot unless none of the
    // facet normals were averaged))

    normals = model.vertexnorms;
    model.vertexnorms = [null];
    for(i = 1 ; i <= model.numvertexnorms ; i++) {
        model.vertexnorms[3 * i + 0] = normals[3 * i + 0];
        model.vertexnorms[3 * i + 1] = normals[3 * i + 1];
        model.vertexnorms[3 * i + 2] = normals[3 * i + 2];
    }

}; // GLvoid

/*
Generates linear texture information given the texture
coordinates in the OBJ file. Gets texture u,v coordinates
*/
objm.loadObj.prototype.objmLinearTexture = function(model) { // (OBJMmodel*)
    var group; // OBJMgroup*
    var dimensions = [0,0,0]; // GLfloat[3]
    var x, y, scalefactor; // GLfloat
    var i; // GLuint

    if(assert(model, 'Model exists')) return;

    // Free texcoords
    // Malloc texcorrds based on numvertices+1

    this.objmDimensions(model,dimensions);
    scalefactor = 2 / this.objmAbs(this.objmMax(this.objmMax(dimensions[0], dimensions[1]), dimensions[2]));

    // Do the calculations
    for(i = i ; i <= model.numvertices ; i++) {
        x = model.vertices[3 * i + 0] * scalefactor;
        y = model.vertices[3 * i + 2] * scalefactor;

        model.texcoords[2 * i + 0] = x / 2;
        model.texcoords[2 * i + 1] = y / 2;
    }

    // Go through and put texture coordinate indices in all the triangles
    group = model.group;
    while(group) {
        for(i = 0 ; i < group.numtriangles() ; i++) {
            objm.T(model,group.triangles[i],'tindices')[0] = objm.T(model,group.triangles[1],'vindices')[0];
            objm.T(model,group.triangles[i],'tindices')[1] = objm.T(model,group.triangles[1],'vindices')[1];
            objm.T(model,group.triangles[i],'tindices')[2] = objm.T(model,group.triangles[1],'vindices')[2];
        }
        group = group.next;
    }
}; // GLvoid

/*
Generates spherical texture information given the texture
coordinates in the OBJ file. Gets spherical coordinates for the texture.
*/
objm.loadObj.prototype.objmSpheremapTexture = function(model) { // (OBJMmodel*)
    var group; // OBJMgroup*
    var theta, phu, rho, x, y, z, r; // GLfloat
    var i; // GLuint

    if(assert(model, 'Model exists')) return;
    if(assert(model.vertexnorms, 'Model has vertex normsl')) return;

    // Free model texcoords if exist
    // Malloc texcoords based on numvertexnorms+1

    for(i = 1 ; i <= model.numvertexnorms ; i++) {
        z = model.vertexnorms[3 * i + 0]; // Re-arrange for pole distortion
        y = model.vertexnorms[3 * i + 1];
        x = model.vertexnorms[3 * i + 2];
        r = Math.sqrt((x *x) + (y * y));
        rho = Math.sqrt((r * r) + (z * z));

        if(r == 0) {
            theta = 0;
            phi = 0;
        } else {
            if(z == 0)
                phi = objm.C_PI / 2;
            else
                phi = Math.acos(z / rho);

            if(y == 0)
                theta = objm.C_PI / 2;
            else
                theta = Math.asin(y / r) + (objm.C_PI / 2);
        }

        model.texcoords[2 * i + 0] = theta / objm.C_PI;
        model.texcoords[2 * i + 1] = phi / objm.C_PI;
    }

    // Go through an put texcoord indices in all the triangles
    group = model.group;
    while(group) {
        for(i = 0 ; i < group.numtriangles ; i++) {
            objm.T(model,group.triangles[i],'tindices')[0] = objm.T(model,group.triangles[i],'nindices')[0];
            objm.T(model,group.triangles[i],'tindices')[1] = objm.T(model,group.triangles[i],'nindices')[1];
            objm.T(model,group.triangles[i],'tindices')[2] = objm.T(model,group.triangles[i],'nindices')[2];
        }
        group = group.next
    }
}; // GLvoid

objm.loadObj.prototype.objmDelete = function(model) { // (OBJMmodel*)
    
}; // GLvoid

/*
Loads a specific OBJ file. Cna be used either
alternatively or together with drawmodel() which
loads a default model and draws it.

TODO: This will have to be updated with ajax server calls.
Models should be saved either in a folder structure or a DB.
*/
objm.loadObj.prototype.objmLoad = function(filename) { // (char*)
    this.pmodel = this.objmReadOBJ(filename);
    if(!this.pmodel) return false;
    this.objmUnitize(this.pmodel);
    this.objmVertexNormals(this.pmodel, 90.0);
    if(!this.pmodel) return false;
    return true;

}; // bool

/*
Saves a loaded OBJ file. The file is saved with
smooth normals, and if material info is present,
a function is called by the writer below to write
out a separate material (.mtl) file. See function
objmWriteOBJ().

TODO: This will have to be updated with ajax server calls.
Models should be saved either in a folder structure or a DB.
*/
objm.loadObj.prototype.objmSave = function(filename) { // (char*)
    if(this.pmodel) {
        this.objmWriteOBJ(this.pmodel, filename, objm.OBJM_SMOOTH_SHADE | objm.OBJM_MATERIAL);
        return true;
    }
    return false
}; // bool

objm.loadObj.prototype.objmReadOBJ = function(filename) { // (char*)
    
}; // OBJMmodel*

objm.loadObj.prototype.objmWriteOBJ = function(model, filename, mode) { // (OBJMmodel*, char*, GLuint)
    
}; // GLvoid

objm.loadObj.prototype.objmDraw = function(model, mode) { // (OBJMmodel*, GLuint)
    
}; // GLvoid

objm.loadObj.prototype.objmList = function(model, mode) { // (OBJMmodel*, GLuint)
    
}; // GLuint

objm.loadObj.prototype.objmWeld = function(model, epsilon) { // (OBJMmodel*, GLfloat)
    
}; // GLvoid

objm.loadObj.prototype.objmReadPPM = function(filename, width, height) { // (char*, int*, int*)
    
}; // GLubyte*

objm.loadObj.prototype.objmMax = function(a, b) { // (GLfloat, GLfloat)
    
}; // GLfloat

objm.loadObj.prototype.objmAbs = function(f) { // (GLfloat)
    
}; // GLfloat

objm.loadObj.prototype.objmDot = function(a, b) { // (GLfloat*, GLfloat*)
    
}; // GLfloat

// Now takes in Vector objects
objm.loadObj.prototype.objmCross = function(a, b, n) { // (GLfloat*, GLfloat*, GLfloat*)
    
}; // GLvoid

objm.loadObj.prototype.objmNormalize = function(v) { // (GLfloat*)
    
}; // GLvoid

objm.loadObj.prototype.objmEqual = function(a, b, epsilon) { // (GLfloat*, GLfloat*, GLfloat)
    
}; // GLboolean

objm.loadObj.prototype.objmWeldVectors = function(vectors, numvectors, epsilon) { // (LGfloat*, GLuint*, GLfloat)
    
}; // GLfloat*

objm.loadObj.prototype.objmFindGroup = function(model, name) { // (OBJMmodel*, char*)
    
}; // OBJMgroup*

objm.loadObj.prototype.objmAddGroup = function(model, name) { // (OBJMmodel*, char*)
    
}; // OBJMgroup*

objm.loadObj.prototype.objmFindMaterial = function(model, name) { // (OBJMmodel*, char*)
    
}; // GLuint

objm.loadObj.prototype.objmDirName = function(path) { // (char*)
    
}; // char*

objm.loadObj.prototype.objmReadMLT = function(model, name) { // (OBJMmodel*, char*)
    
}; // GLvoid

objm.loadObj.prototype.objmWriteMLT = function(model, modelpath, mtllibname) { // (OBJMmodel*, char*, char*)
    
}; // GLvoid

objm.loadObj.prototype.objmFirstPass = function(model, file) { // (OBJMmodel*, FILE*)
    
}; // GLvoid

objm.loadObj.prototype.objmSecondPass = function(model, file) { // (OBJMmodel*, FILE*)
    
}; // GLvoid

objm.loadObj.prototype.normalize = function(v) { // (float*)
    
}; // float

objm.loadObj.prototype.invert = function(src, inverse) { // (GLdouble[16], GLdouble[16])
    
}; // GLboolean

objm.loadObj.prototype.identity = function(m) { // (GLdouble[16])
    
}; // void

/*
Loads a default OBJ file. The path to the file may be changed
as desired. Loads up material from the approppriate material file if
so specified in the OBJ file. See function objmDraw().
*/
objm.loadObj.prototype.drawModel = function() { // ()
    if(!this.pmodel) {
        this.pmodel = this.objmReadOBJ("default.obj");
        if(!this.pmodel) {
            console.error("Default model 'default.obj' not found by objmReadOBJ");
            return;
        }
        this.objmUnitize(this.pmodel);
        this.objmFacetNormals(this.pmodel);
        this.objmVertexNormals(this.pmodel, 90.0);
    }

    this.objmDraw(this.pmodel, objm.OBJM_SMOOTH_SHADE | objm.OBJM_MATERIAL);
}; // void

objm.loadObj.prototype.buildEdgeList = function(model, numtriangles) { // (OBJMmodel*, int)
    
}; // void



/*********************************************
**********************************************/

/*
Assert Function - Found on NetTuts and modified
http://net.tutsplus.com/tutorials/javascript-ajax/quick-tip-quick-and-easy-javascript-testing-with-assert/

If there is an assert, meaning outcome is false, null or 0, print an error
message to the console and return true. If there is no assert, meaning
outcome has a valid value, return false, as no assert occured.
*/

function assert(outcome, description) {
    if(!outcome) console.error("Failed Assert \""+description+"\"");
    return !outcome;
}