/*!
 * Bootstrap Context Menu
 * Version: 2.2
 * A small variation of the dropdown plugin by @Festum
 * https://github.com/Festum/bootstrap-contextmenu
 *
 * New options added by @jeremyhubble for javascript launching
 *  $('#elem').contextmenu({target:'#menu',before:function(e) { return true; } });
 *
 *
 * Bootstrap (http://getbootstrap.com/).
 */

/* =========================================================
 * bootstrap-contextmenu.js
 * =========================================================
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ========================================================= */
(function ($) {
	var j = {};
	var k = {
		init: function (a, b) {
			var c;
			c = k.guid();
			var d = $('<div class="dropdown clearfix bootstrap-contextmenu"></div>');
			d.css({
				display: 'none',
				position: 'absolute'
			}).attr({
				id: c
			}).appendTo('body');
			$(a).data('e_id', c);
			mc = $(a).data('mc');
			if (mc == undefined) {
				mc = 3
			}
			$(a).data('e', $.extend(true, {}, b));
			$(a).css('-moz-user-select', 'none');
			$(a).css('-khtml-user-select', 'none');
			$(a).css('user-select', 'none');
			k.buildElements(a, d, b, true);
			k.create(a, mc)
		},
		buildElements: function (a, b, c, d) {
			var f = $('<ul aria-labelledby="dropdownMenu" role="menu" class="dropdown-menu"></ul>');
			d = (d != undefined) ? d : true;
			if (d == true) {
				f.show()
			}
			for (var i in c) {
				var g = $('<li></li>');
				g.attr({
					id: i
				});
				switch (typeof c[i]) {
				case 'object':
					if (c[i].text != undefined) {
						if (c[i].text == '---') {
							g.addClass('divider')
						} else {
							var h = $('<a href="#" tabindex="-1"></a>');
							if (typeof c[i].icon == 'string') {
								g.append(h.append(' <i class="' + c[i].icon + '"></i> ').append(c[i].text))
							} else {
								h = $('<a href="#" tabindex="-1">' + c[i].text + '</a>');
								g.append(h)
							} if (typeof c[i].click == 'function') {
								h.bind('mousedown', {
									key: i,
									target: a,
									callback: c[i].click
								}, function (e) {
									if (!$(this).parent().hasClass('disabled')) {
										e.data.callback(e.data.target, $(this).parent())
									}
								})
							}
							if (c[i].disabled != undefined && c[i].disabled == true) {
								g.addClass('disabled')
							}
							if (typeof c[i].children == 'object') {
								if (!g.hasClass('disabled')) {
									g.addClass('dropdown-submenu');
									k.buildElements(a, g, c[i].children, false)
								}
							}
							if (g.hasClass('disabled')) {
								g.children('a').children('i').hide()
							}
						}
					}
					break;
				case 'string':
					if (c[i] == '---') {
						g.addClass('divider')
					} else {
						g.append($('<a href="#" tabindex="-1">' + c[i] + '</a>'))
					}
					break
				}
				f.append(g)
			}
			b.append(f)
		},
		create: function (f, g) {
			if (g == undefined) g = 3;
			if (/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent)) {
				g = 1
			}
			var h = $(f).data('e_id');
			$("#" + h).bind("contextmenu", function (e) {
				e.preventDefault()
			});
			$(f).mousedown(function (e) {
				if (!(/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent))) {
					$("#" + h).hide()
				}
			});
			$(f).bind("contextmenu", function (e) {
				e.preventDefault()
			});
			$(document).mousedown(function (e) {
				if (!(/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent))) {
					$("#" + h).hide()
				}
			});
			if (/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent)) {
				$("#" + h).mousedown(function (e) {
					$("#" + h).hide()
				});
				$(document).mousedown(function (e) {
					$('.bootstrap-contextmenu').hide()
				})
			}
			$(f).mousedown(function (e) {
				if (!(/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent))) {
					$('.bootstrap-contextmenu').hide()
				} else {
					$('.bootstrap-contextmenu').each(function () {
						if ($(this).attr('id') != h) {
							$(this).hide()
						}
					})
				}
				var c = e;
				e.stopPropagation();
				$(this).mouseup(function (e) {
					e.stopPropagation();
					$(this).unbind('mouseup');
					if (c.which == g) {
						if (j[h] != undefined) {
							var a = j[h];
							a($(f))
						}
						var d = {}, x, y;
						if (self.innerHeight) {
							d.pageYOffset = self.pageYOffset;
							d.pageXOffset = self.pageXOffset;
							d.innerHeight = self.innerHeight;
							d.innerWidth = self.innerWidth
						} else if (document.documentElement && document.documentElement.clientHeight) {
							d.pageYOffset = document.documentElement.scrollTop;
							d.pageXOffset = document.documentElement.scrollLeft;
							d.innerHeight = document.documentElement.clientHeight;
							d.innerWidth = document.documentElement.clientWidth
						} else if (document.body) {
							d.pageYOffset = document.body.scrollTop;
							d.pageXOffset = document.body.scrollLeft;
							d.innerHeight = document.body.clientHeight;
							d.innerWidth = document.body.clientWidth
						}(e.pageX) ? x = e.pageX : x = e.clientX;
						(e.pageY) ? y = e.pageY : y = e.clientY;
						var b = $("#x-context-" + $(f).attr('id')).height();
						if (/Android|webOS|iPhone|iPad|iPod|BlackBerry/i.test(navigator.userAgent)) {
							if ($("#" + h).css('display') == 'block') {
								$('.bootstrap-contextmenu').hide()
							} else {
								if (y + b > $(document).height()) {
									$("#" + h).css({
										top: y - b,
										left: x
									}).fadeIn(20)
								} else {
									$("#" + h).css({
										top: y,
										left: x
									}).fadeIn(20)
								}
							}
						} else {
							if (y + b > $(document).height()) {
								$("#" + h).css({
									top: y - b,
									left: x
								}).fadeIn(20)
							} else {
								$("#" + h).css({
									top: y,
									left: x
								}).fadeIn(20)
							}
						}
					}
				})
			})
		},
		guid: function () {
			var a = function () {
				return Math.floor(Math.random() * 0x10000).toString(16)
			};
			return (a() + a() + "-" + a() + "-" + a() + "-" + a() + "-" + a() + a() + a())
		}
	};
	$.fn.contextMenu = function (a, b) {
		b = (b != undefined) ? b : '';
		if (typeof a == 'string' && b != '') {
			var c = $(this).data('e_id');
			var d = $(this).data('e');
			if (typeof c == 'string' && c != '') {
				if (typeof b != 'function') {
					var e = [];
					if (b.indexOf('>') > 0) {
						var f = b.split('>');
						for (var i in f) {
							e[e.length] = $.trim(f[i])
						}
					} else {
						e[e.length] = $.trim(b)
					}
				}
				switch (a) {
				case 'disable':
					$('#' + c).remove();
					var g = d;
					var l = e.length - 1;
					for (var i in e) {
						if (g[e[i]] != undefined) {
							if (i == l) {
								g[e[i]].disabled = true
							} else {
								if (g[e[i]].children != undefined) {
									g = g[e[i]].children
								} else {
									break
								}
							}
						} else {
							break
						}
					}
					k.init(this, d);
					break;
				case 'enable':
					$('#' + c).remove();
					var g = d;
					var l = e.length - 1;
					for (var i in e) {
						if (g[e[i]] != undefined) {
							if (i == l) {
								g[e[i]].disabled = false
							} else {
								if (g[e[i]].children != undefined) {
									g = g[e[i]].children
								} else {
									break
								}
							}
						} else {
							break
						}
					}
					k.init(this, d);
					break;
				case 'beforeDisplay':
					if (typeof b == 'function') {
						j[c] = b
					}
					break
				}
			}
		} else if (typeof a == 'object') {
			return $(this).each(function () {
				k.init(this, a)
			})
		}
	}
})(jQuery);