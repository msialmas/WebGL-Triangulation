var modelViewMatrix;
var projectionMatrix;
var modelViewMatrixStack;

function createGLContext(canvas) {
    var names = ["webgl", "experimental-webgl"];
    var context = null;
    for (var i = 0; i < names.length; i++) {
        try {
            context = canvas.getContext(names[i]);
        } catch (e) {
        }
        if (context) {
            break;
        }
    }
    if (context) {
        context.viewportWidth = canvas.width;
        context.viewportHeight = canvas.height;
    } else {
        alert("Failed to create WebGL context!");
    }
    return context;
}

function loadShaderFromDOM(id) {
    var shaderScript = document.getElementById(id);

    // If we don't find an element with the specified id
    // we do an early exit 
    if (!shaderScript) {
        return null;
    }

    // Loop through the children for the found DOM element and
    // build up the shader source code as a string
    var shaderSource = "";
    var currentChild = shaderScript.firstChild;
    while (currentChild) {
        if (currentChild.nodeType == 3) { // 3 corresponds to TEXT_NODE
            shaderSource += currentChild.textContent;
        }
        currentChild = currentChild.nextSibling;
    }

    var shader;
    if (shaderScript.type == "x-shader/x-fragment") {
        shader = gl.createShader(gl.FRAGMENT_SHADER);
    } else if (shaderScript.type == "x-shader/x-vertex") {
        shader = gl.createShader(gl.VERTEX_SHADER);
    } else {
        return null;
    }

    gl.shaderSource(shader, shaderSource);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert(gl.getShaderInfoLog(shader));
        return null;
    }
    return shader;
}

function setupShaders() {
    var vertexShader = loadShaderFromDOM("shader-vs");
    var fragmentShader = loadShaderFromDOM("shader-fs");

    shaderProgram = gl.createProgram();
    gl.attachShader(shaderProgram, vertexShader);
    gl.attachShader(shaderProgram, fragmentShader);
    gl.linkProgram(shaderProgram);

    if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
        alert("Failed to setup shaders");
    }

    gl.useProgram(shaderProgram);
    var resolutionLocation = gl.getUniformLocation(shaderProgram, "u_resolution");
    gl.uniform2f(resolutionLocation, canvas.width, canvas.height);
    var zValue = gl.getUniformLocation(shaderProgram, "u_ZValue");
    gl.uniform1f(zValue, canvas.width);

    u_NormalMatrix = gl.getUniformLocation(shaderProgram, 'u_NormalMatrix');
    u_LightColor = gl.getUniformLocation(shaderProgram, 'u_LightColor');
    u_LightPosition = gl.getAttribLocation(shaderProgram, "a_LightPosition");
    u_AmbientLight = gl.getUniformLocation(shaderProgram, 'u_AmbientLight');
}

function pushModelViewMatrix() {
    var copyToPush = new Matrix4(g_MvpMatrix);
    modelViewMatrixStack.push(copyToPush);
}

function popModelViewMatrix() {
    if (modelViewMatrixStack.length == 0) {
        throw "Error popModelViewMatrix() - Stack was empty ";
    }
    g_MvpMatrix = modelViewMatrixStack.pop();
}


function uploadModelViewMatrixToShader() {
    gl.uniformMatrix4fv(shaderProgram.uniformMVMatrix, false, modelViewMatrix);
}

function uploadProjectionMatrixToShader() {
    gl.uniformMatrix4fv(shaderProgram.uniformProjMatrix, false, projectionMatrix);
}

function hex2num(hex) {
    if (hex.charAt(0) == "#")
        hex = hex.slice(1); //Remove the '#' char - if there is one.
    hex = hex.toUpperCase();
    var hex_alphabets = "0123456789ABCDEF";
    var value = new Array(3);
    var k = 0;
    var int1, int2;
    for (var i = 0; i < 6; i += 2) {
        int1 = hex_alphabets.indexOf(hex.charAt(i));
        int2 = hex_alphabets.indexOf(hex.charAt(i + 1));
        value[k] = (int1 * 16) + int2;
        k++;
    }
    return(value);
}
function initTexture(filename) {
    // Create a texture object
    var texture = gl.createTexture();
    if (!texture) {
        console.log('Failed to create the texture object');
        return null;
    }

    // Get the storage location of u_Sampler
    var u_Sampler = gl.getUniformLocation(shaderProgram, 'u_Sampler');
    if (!u_Sampler) {
        console.log('Failed to get the storage location of u_Sampler');
        return null;
    }

    // Create the image object
    var image = new Image();
    if (!image) {
        console.log('Failed to create the image object');
        return null;
    }
    // Register the event handler to be called when image loading is completed
    image.onload = function() {
        loadTexture(texture, u_Sampler, image);
    };
    // Tell the browser to load an Image
    image.src = 'resources/' + filename;

    return texture;
}
function loadTexture(texture, u_Sampler, image) {
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, 1);  // Flip the image Y coordinate
    // Activate texture unit0
    gl.activeTexture(gl.TEXTURE0);
    // Bind the texture object to the target
    gl.bindTexture(gl.TEXTURE_2D, texture);

    // Set texture parameters
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    // Set the image to texture
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB, gl.RGB, gl.UNSIGNED_BYTE, image);

    // Pass the texure unit 0 to u_Sampler
    gl.uniform1i(u_Sampler, 0);
}
function initEventHandlers(canvas, currentAngle) {
    var dragging = false;         // Dragging or not
    var lastX = -1, lastY = -1;   // Last position of the mouse

    canvas.onmousedown = function(ev) {   // Mouse is pressed
        var x = ev.clientX, y = ev.clientY;
        // Start dragging if a moue is in <canvas>
        var rect = ev.target.getBoundingClientRect();
        if (rect.left <= x && x < rect.right && rect.top <= y && y < rect.bottom) {
            lastX = x;
            lastY = y;
            dragging = true;
        }
    };

    canvas.onmouseup = function(ev) {
        dragging = false;
    }; // Mouse is released

    canvas.onmousemove = function(ev) { // Mouse is moved
        var x = ev.clientX, y = ev.clientY;
        if (dragging) {
            var factor = 100 / canvas.height; // The rotation ratio
            var dx = factor * (x - lastX);
            var dy = factor * (y - lastY);
            // Limit x-axis rotation angle to -90 to 90 degrees
            currentAngle[0] = Math.max(Math.min(currentAngle[0] + dy, 90.0), -90.0);
            currentAngle[1] = currentAngle[1] + dx;
        }
        lastX = x, lastY = y;
    };
}


function handleFileSelect(evt) {
    var files = evt.target.files; // FileList object
    return files;
}


function getTexture(schema) {
    var file = getFileInfo();
    schema.textureData = initTexture(file);
    schema.isColorMode = false;
}

function getFileInfo() {
    var fileInput = document.getElementById("fileInput");
    if ('files' in fileInput) {
        if (fileInput.files.length > 0) {
            for (var i = 0; i < fileInput.files.length; i++) {
                return fileInput.files[i].name;
            }
        }
    }
}



function getColor()
{
    var color = hex2num(document.getElementById("colorpicker").value);
    return getRGB(color);
}
function getLightColor()
{
    var color = hex2num(document.getElementById("lightcolorpicker").value);
    return getRGB(color);
}
function getAmbientColor()
{
    var color = hex2num(document.getElementById("ambientcolorpicker").value);
    return getRGB(color);
}

function getLightPosition() {
    var x = document.getElementById("lightx").value;
    var y = document.getElementById("lighty").value;
    var z = document.getElementById("lightz").value;
    var xyz = new Array([x], [y], [z]);
    return xyz;
}

function getRGB(color) {
    var binR = color[0] / 255;
    var binG = color[1] / 255;
    var binB = color[2] / 255;
    var rgb = new Array([binR], [binG], [binB]);
    return rgb;
}

function numToHex(num) {
    return ("00" + num.toString(16).toUpperCase()).slice(-2);
}
function setColorPicker(r, g, b) {
    var decR = r * 255;
    var decG = g * 255;
    var decB = b * 255;
    var hexR = numToHex(decR);
    var hexG = numToHex(decG);
    var hexB = numToHex(decB);
    document.getElementById("colorpicker").value = "#" + hexR + hexG + hexB;
}

var modify = {
    COLOR: 0,
    TEXTRURE: 1,
    SIZE: 2,
    NEW: 3
};