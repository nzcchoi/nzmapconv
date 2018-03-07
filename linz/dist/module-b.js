"use strict";
var linz;
(function (linz) {
    var moduleB;
    (function (moduleB) {
        var B1 = /** @class */ (function () {
            function B1() {
            }
            B1.prototype.test = function () {
                console.log('B1-test');
            };
            return B1;
        }());
        moduleB.B1 = B1;
    })(moduleB = linz.moduleB || (linz.moduleB = {}));
})(linz || (linz = {}));
//# sourceMappingURL=module-b.js.map