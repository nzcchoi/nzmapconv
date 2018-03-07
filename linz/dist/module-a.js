"use strict";
var linz;
(function (linz) {
    var moduleA;
    (function (moduleA) {
        var A1 = /** @class */ (function () {
            function A1() {
            }
            A1.prototype.test = function () {
                console.log('A1-test');
            };
            return A1;
        }());
        moduleA.A1 = A1;
    })(moduleA = linz.moduleA || (linz.moduleA = {}));
})(linz || (linz = {}));
//# sourceMappingURL=module-a.js.map