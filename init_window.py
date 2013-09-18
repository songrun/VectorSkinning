#!/opt/local/bin/python
from Tkinter import *  
import ttk as ttk
import tkMessageBox
import numpy as np
from weight_inverse_distance import *
from bezier_chain_constraints import *

gOnce = False

class Window:
	canvas = None
	mode = None #indicator of one of the three modes.
	constraint = None #indicator of which constraint is being used.
	root = None
	selected = None
	popup = None
	traceT, traceR, traceS = [], [], []
	transforms = {}
	
	def __init__(self, parent):
	
		self.root = parent
		mainWindow = ttk.Frame(parent)	
		self.init_UI(mainWindow)	
		mainWindow.grid()
		return
		
	def init_UI(self, parent):
	
		menubar = ttk.Frame(parent, relief=GROOVE, borderwidth=1, width=800, height=28)
		self.init_menubar(menubar)
		
		self.canvas = Canvas(parent, width=800, height=600, bd=2, cursor='pencil', relief=SUNKEN)
		self.canvas.bind("<Button-1>", self.onclick_handler)
		self.canvas.bind("<Shift-B1-Motion>", self.on_shift_mouse_handler)
		self.canvas.bind("<Control-B1-Motion>", self.on_control_mouse_handler)
#		self.canvas.bind("<Double-Button-1>", self.on_double_click_handler)
		self.canvas.bind("<B1-Motion>", self.on_motion_handler)
		self.canvas.bind("<ButtonRelease-1>", self.onrelease_handler)
		self.canvas.grid()
		return	
	'''
	A series of actions a user can do
	'''		
	def onclick_handler(self, event):
	
		mode = self.mode.get()
		if mode == 0:
			controls = self.canvas.find_withtag('controls')
			if len( controls ) == 4:
				self.canvas.delete( controls[0] )
			
			r= 3
			x0, x1, y0, y1 = event.x-r, event.x+r, event.y-r, event.y+r
			self.canvas.create_oval(x0, y0, x1, y1, fill='blue', outline='blue', tags='controls')
			
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_bezier_curve( )
				if len(self.canvas.find_withtag('handles')) > 0:
					self.redraw_handle_affected_curve()
					self.redraw_approximated_bezier()

		elif mode == 1:
			r = 3
			x0, x1, y0, y1 = event.x-r, event.x+r, event.y-r, event.y+r
			handle = self.canvas.create_rectangle(x0, y0, x1, y1, fill='red', outline='red', tags='handles')
			self.transforms[handle] = np.identity(3).reshape(-1)
			self.popup_handle_editor( handle )
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()

		elif mode == 2:
			overlaps = self.canvas.find_overlapping(event.x-3, event.y-3, event.x+3, event.y+3)
			sels = set(overlaps) & set(self.canvas.find_withtag('handles')) 
			if len(sels) > 0:
				self.selected = sels.pop()
				self.popup_handle_editor(self.selected)

	# for translation ---- press and drag
	def on_motion_handler(self, event):
	
		if self.selected != None and self.popup != None: 
			h = self.selected
			if len(self.traceT) != 2:
				self.traceT = [event.x, event.y]
				return
				
			trans_x, trans_y = event.x-self.traceT[0], event.y-self.traceT[1]
			self.traceT = [event.x, event.y]
			
			window_name = self.canvas.itemcget( self.popup, 'window' )
			p = self.root.nametowidget( window_name )

			labelFrame = p.winfo_children()[0]
			entry_transX = labelFrame.winfo_children()[2]
			entry_transY = labelFrame.winfo_children()[5]
			
			ori_x, ori_y = self.transforms[h][2], self.transforms[h][5]
			new_x, new_y = trans_x+ori_x, trans_y+ori_y
			entry_transX.delete(0,END)
			entry_transX.insert(0,new_x)
			entry_transY.delete(0,END)
			entry_transY.insert(0,new_y)
			self.transforms[h][2], self.transforms[h][5] = new_x, new_y
			
#			self.canvas.move(self.selected, event.x-ori_x, event.y-ori_y)
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()				
		return
	
	# for rotation	---- Shift + press and drag
	def on_shift_mouse_handler(self, event):
	
		if self.selected != None and self.popup != None: 
			h = self.selected
			if len(self.traceR) != 2:
				coord = self.canvas.coords( h )
				self.traceR = [(coord[0]+coord[2])/2, (coord[1]+coord[3])/2]
				
			coord = self.canvas.coords( h )
			origin = [(coord[0]+coord[2])/2, (coord[1]+coord[3])/2]
			vec_0 = array([self.traceR[0]-origin[0], self.traceR[1]-origin[1]])
			vec_1 = array([event.x-origin[0], event.y-origin[1]])
			self.traceR = [event.x, event.y]
			
			assert len(vec_0) == len(vec_1)
			len_0 = sum([vec_0[i]**2 for i in range(len(vec_0))])**.5
			len_1 = sum([vec_1[i]**2 for i in range(len(vec_1))])**.5
			if len_0*len_1 == 0:
				return
			costheta = dot(vec_1, vec_0)/(len_0*len_1)
			sintheta = cross(vec_1, vec_0)/(len_0*len_1)
			
			window_name = self.canvas.itemcget( self.popup, 'window' )
			p = self.root.nametowidget( window_name )

			labelFrame = p.winfo_children()[0]
			entry_1 = labelFrame.winfo_children()[0]
			entry_2 = labelFrame.winfo_children()[1]
			entry_3 = labelFrame.winfo_children()[3]
			entry_4 = labelFrame.winfo_children()[4]
			
			entry_5 = labelFrame.winfo_children()[2]
			entry_6 = labelFrame.winfo_children()[5]
			
			R = array([[costheta, sintheta, 0], [-sintheta, costheta, 0], [0, 0, 1]])
			T = array([[1, 0, -origin[0]], [0, 1, -origin[1]], [0, 0, 1]])

			newM = dot((dot( dot(np.linalg.inv(T), R), T )), self.transforms[h].reshape(3, -1))
			self.transforms[h] = newM.reshape(-1)			
							
			entry_1.delete(0,END)
			entry_1.insert(0,newM[0][0])
			entry_2.delete(0,END)
			entry_2.insert(0,newM[0][1])
			entry_3.delete(0,END)
			entry_3.insert(0,newM[1][0])
			entry_4.delete(0,END)
			entry_4.insert(0,newM[1][1])
			entry_5.delete(0,END)
			entry_5.insert(0,newM[0][2])
			entry_6.delete(0,END)
			entry_6.insert(0,newM[1][2])
			
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()				
		return
	
	# for scaling ---- Control + press and drag
	def on_control_mouse_handler(self, event):
	
		if self.selected != None and self.popup != None:
			h = self.selected 
			if len(self.traceS) != 2:
				self.traceS = [event.x, event.y]
				return
				
			coord = self.canvas.coords( h )
			origin = [(coord[0]+coord[2])/2, (coord[1]+coord[3])/2]

			box = self.canvas.bbox('original_bezier')
			width, height = box[2]-box[0], box[3]-box[1]
			scaleX, scaleY = np.power(2., float(event.x-self.traceS[0])/width), np.power(2., float(event.y-self.traceS[1])/height)
			
			self.traceS = [event.x, event.y]
			
			T = array([[1, 0, -origin[0]], [0, 1, -origin[1]], [0, 0, 1]])
			S = array([[scaleX, 0., 0.], [0., scaleY, 0.], [0., 0., 1.]])
 			newM = dot( np.linalg.inv(T), dot( S, dot( T, self.transforms[h].reshape(3, -1) )) )

			self.transforms[h] = newM.reshape(-1)

			
			window_name = self.canvas.itemcget( self.popup, 'window' )
			p = self.root.nametowidget( window_name )

			labelFrame = p.winfo_children()[0]
			entry_1 = labelFrame.winfo_children()[0]
			entry_2 = labelFrame.winfo_children()[1]
			entry_3 = labelFrame.winfo_children()[3]
			entry_4 = labelFrame.winfo_children()[4]
			
			entry_5 = labelFrame.winfo_children()[2]
			entry_6 = labelFrame.winfo_children()[5]
			entry_1.delete(0,END)
			entry_1.insert(0,newM[0][0])
			entry_2.delete(0,END)
			entry_2.insert(0,newM[0][1])
			entry_3.delete(0,END)
			entry_3.insert(0,newM[1][0])
			entry_4.delete(0,END)
			entry_4.insert(0,newM[1][1])
			entry_5.delete(0,END)
			entry_5.insert(0,newM[0][2])
			entry_6.delete(0,END)
			entry_6.insert(0,newM[1][2])
			
			if len( self.canvas.find_withtag('controls') ) == 4:
				self.redraw_handle_affected_curve()
				self.redraw_approximated_bezier()		
						
		return	
	
	def onrelease_handler(self, event):
	
		self.traceT = []
		self.traceR = []
		self.traceS = []
		
		if self.selected == None:
			self.canvas.delete('popup') 
		else:	
			self.selected = None

	'''
		Actions End
	'''
	'''
		algorithms for drawing
	'''		
	def redraw_bezier_curve(self):
	
		self.canvas.delete( 'original_bezier' )
		cps = self.canvas.find_withtag('controls')
		cps = [self.canvas.bbox(x) for x in cps]
		cps = [((x[0]+x[2])/2, (x[1]+x[3])/2) for x in cps]
		cps = np.mat(cps)
		known = np.asarray( M*cps ).reshape(4, -1)
		
		ps = [] 
		for t in range(0, 101):
			p = np.dot( np.asarray( [(float(t)/100)**3, (float(t)/100)**2, float(t)/100, 1] ), known )
			ps = ps + [p[0], p[1]]	
		self.canvas.create_line(ps, smooth=True, width=2, tags='original_bezier')
		
	def redraw_handle_affected_curve(self):
	
		self.canvas.delete( 'handle_affected' )
		pts = self.canvas.find_withtag( 'original_bezier' )
		pts = self.canvas.coords( pts )
		
		tps = []
		for i in range(0, len(pts)/2):
			p = np.asarray([pts[i*2], pts[i*2+1], 1.])
			m = np.zeros(9)
			w = {}
			for h in self.canvas.find_withtag('handles'):
			# coefficient decided by square of distance
				pos = self.canvas.coords(h)
				w[h] = 1 / ( (pos[0]/2+pos[2]/2-p[0])**2 + (pos[1]/2+pos[3]/2-p[1])**2 )
				
			for h in self.canvas.find_withtag('handles'):
				m = m + self.transforms[h]*( w[h] / sum(w.values()) )
			
			p = np.dot( m.reshape(3, 3), p.reshape(3,-1) ).reshape(-1)
			tps = tps + [p[0], p[1]]

		self.canvas.create_line(tps, smooth=True, width=2, fill='magenta', tags='handle_affected')
		
	def redraw_approximated_bezier(self):
	
		self.canvas.delete( 'approximated' )
		self.canvas.delete( 'new_controls' )
		
		handles = []
		trans = []
		H_set = self.canvas.find_withtag('handles')
		for i in range( len( H_set ) ):
			pos = self.canvas.coords( H_set[i] )
			pos = [(pos[0]+pos[2])/2, (pos[1]+pos[3])/2, 1.]
			handles.append(pos)
			trans.append( self.transforms[H_set[i]] )
		
		cps = self.canvas.find_withtag('controls')
		cps = [self.canvas.bbox(x) for x in cps]
		cps = [((x[0]+x[2])/2, (x[1]+x[3])/2, 1) for x in cps]
		cps = np.mat(cps)
		
		'''
		global gOnce
		
		import csv
		outname = 'integration_accuracy.csv'
		def uniquepath( path ):
			import os
			i = 1
			result = path
			while os.path.exists( result ):
				split = os.path.splitext( path )
				result = split[0] + ( ' %d' % i ) + split[1]
				i += 1
			return result
		
		writer = csv.writer( open( uniquepath( 'integration_accuracy.csv' ), 'w' ) )
		'''
		
# 		partition = [0., 0.2, 0.5, 1.]
		partition = [0.2, 0.8]
# 		partition = [0., 0.1, 1.]

		assert sum( partition ) == 1.0
		for x in partition:
			assert x > 0. and x <=1.
		
		level = self.constraint.get()
		P_primes = compute_control_points_chain_with_constraint( partition, cps, handles, trans, level)

		#del writer
		colors = ['green', 'gold', 'brown', 'coral', 'cyan', 'gray']
		# new control points
		for i in range( len( P_primes ) ):
			for pp in P_primes[i]:
				pp = np.asarray( pp ).reshape(3)
				r= 3
				x0, x1, y0, y1 = pp[0]-r, pp[0]+r, pp[1]-r, pp[1]+r
				self.canvas.create_oval(x0, y0, x1, y1, fill=colors[i], outline=colors[i], tags='new_controls')

		# new bezier curve
# 		num_sample = 100
# 		for i in range( len( P_primes ) ):
# 			known = np.asarray( M * P_primes[i] ).reshape(4, -1)	
# 			ps = [] 
# 			samples = array( range( num_sample ) ) / float(num_sample - 1)
# 			
# 			start, end = np.where( samples >= sections[i] )[0][0], np.where( samples <= sections[i+1] )[0][-1]
# 			print 'start , end ', start, end, 'hehe', np.where( samples >= sections[i] )
# 			for t in samples[start: end]:
# 				p = np.dot( np.asarray( [t**3, t**2, t, 1] ), known )
# 				ps = ps + [p[0], p[1]]	
# 			self.canvas.create_line(ps, smooth=True, width=2, fill=colors[i]+'1', tags='approximated')	 
		num_sample = 100
		for i in range( len( P_primes ) ):
			known = np.asarray( M * P_primes[i] ).reshape(4, -1)	
			ps = [] 
			samples = array( range( num_sample ) ) / float(num_sample - 1)
			
# 			start, end = np.where( samples >= sections[i] )[0][0], np.where( samples <= sections[i+1] )[0][-1]
			for t in samples:
				p = np.dot( np.asarray( [t**3, t**2, t, 1] ), known )
				ps = ps + [p[0], p[1]]	
			self.canvas.create_line(ps, smooth=True, width=2, fill=colors[i]+'1', tags='approximated')	 

	'''
		Drawing End
	'''	
	def init_menubar(self, menubar):
	
		menubar.grid_propagate(0)
		menubar.grid()
		
		#A menu in Tk is a combination of a Menubutton (the title of the
		#menu) and the Menu (what drops down when the Menubutton is pressed)
			
		mb_file = ttk.Menubutton(menubar, text='file')
		mb_file.grid(column=0, row=0, ipadx=20)
		mb_file.menu = Menu(mb_file)
		
		#Once we've specified the menubutton and the menu, we can add
		#different commands to the menu
		
		mb_file.menu.add_command(label='open')
		mb_file.menu.add_command(label='close')
		
		mb_edit = ttk.Menubutton(menubar, text='edit')
		mb_edit.grid(column=1, row=0, ipadx=20 )
		mb_edit.menu = Menu(mb_edit)
		
		self.mode = IntVar()
		mb_edit.menu.add_radiobutton(label='add control point', variable=self.mode, value=0, command=self.change_mode)
		mb_edit.menu.add_radiobutton(label='add handle', variable=self.mode, value=1, command=self.change_mode)
		mb_edit.menu.add_radiobutton(label='edit handle', variable=self.mode, value=2, command=self.change_mode)
		mb_edit.menu.add_separator()
		mb_edit.menu.add_command(label='delete controls', command=self.del_controls)
		mb_edit.menu.add_command(label='delete handles', command=self.del_handles)
		mb_edit.menu.add_command(label='clear canvas', command=self.clear_canvas)
		
		
		## menu dealing with constrains 
		mb_cons = ttk.Menubutton(menubar, text='constrains')
		mb_cons.grid(column=2, row=0)
		mb_cons.menu = Menu(mb_cons)
		self.constraint = IntVar()
		mb_cons.menu.add_radiobutton(label='No constrains', variable=self.constraint, value=0, command=self.change_constrain)		
		mb_cons.menu.add_radiobutton(label='C^0', variable=self.constraint, value=1, command=self.change_constrain)
		mb_cons.menu.add_radiobutton(label='C^1', variable=self.constraint, value=2, command=self.change_constrain)
		mb_cons.menu.add_radiobutton(label='G^1', variable=self.constraint, value=3, command=self.change_constrain)
		
		mb_help = ttk.Menubutton(menubar,text='help')
		mb_help.grid(column=3, row=0, padx=300, sticky=E)
		
		mb_file['menu'] = mb_file.menu
		mb_edit['menu'] = mb_edit.menu
		mb_cons['menu'] = mb_cons.menu
# 		mb_help['menu'] = mb_help.menu
		return 
		
	def popup_handle_editor(self, handle_id):
	
		self.canvas.delete('popup') 
		self.selected = handle_id
		
		coord = self.canvas.bbox(handle_id)
		data = self.transforms[handle_id]
		values = []
		for i in range(0, 9):
			values.append(StringVar())
			values[i].set(data[i])
		
		frame = Frame(self.root, borderwidth=1, relief=RIDGE)
		labelFrame = LabelFrame(frame, text='Transform', relief=GROOVE, borderwidth=1)
		labelFrame.grid()
		w11 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[0])
		w11.grid(row=0, column=0)
		w12 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[1])
		w12.grid(row=0, column=1)
		w13 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[2])
		w13.grid(row=0, column=2)
		w21 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[3])
		w21.grid(row=1, column=0)
		w22 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[4])
		w22.grid(row=1, column=1)
		w23 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[5])
		w23.grid(row=1, column=2)
		w31 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[6])
		w31.grid(row=2, column=0)
		w32 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[7])
		w32.grid(row=2, column=1)
		w33 = Entry(labelFrame, cursor='xterm', highlightcolor='cyan', width=6, textvariable=values[8])
		w33.grid(row=2, column=2)
		
		popup = self.canvas.create_window((coord[0]+coord[2])/2, (coord[1]+coord[3])/2, anchor=NW, width=204, height=140, window=frame, tags='popup')	
		
		btn1 = Button(frame, text='Save', command=lambda id=handle_id, vals=values, popup=popup : self.save_transforms(id, vals, popup), width=4)
		btn1.grid(row=3, column=0, sticky=SW)
		btn2 = Button(frame, text='Delete', command=lambda id=handle_id, popup=popup : self.delete_handle(id, popup), width=4)
		btn2.grid(row=3, column=0, sticky=S)
		btn3 = Button(frame, text='Cancel', command=lambda popup=popup : self.remove_popup(popup), width=4)
		btn3.grid(row=3, column=0, sticky=SE)
	
		self.popup = popup
#		return popup



	def onExit(self):
	
		self.quit()


	def change_mode(self):
	
		mode = self.mode.get()
		
		if mode == 0:
			self.canvas.config(cursor = 'pencil')
		elif mode == 1: 
			self.canvas.config(cursor = 'target')
		elif mode == 2:
			self.canvas.config(cursor = 'hand1')
			
	def change_constrain(self):
		
		if ( len( self.canvas.find_withtag( 'approximated' ) ) != 0 ):
			self.redraw_approximated_bezier()

	def del_controls(self):
	
		self.canvas.delete('controls')
		self.canvas.delete('original_bezier')
		self.canvas.delete('handle_affected')
		self.canvas.delete('approximated')
		
	
	def del_handles(self):
	
		self.canvas.delete('handles')
		self.canvas.delete('handle_affected')
		self.canvas.delete('approximated')
		self.transforms.clear()

	def clear_canvas(self):
	
		self.canvas.delete('all')
		self.transforms.clear()
	
	def save_transforms(self, id, vals, popup):
		
		for i in range(0, 9):
			self.transforms[id][i] = vals[i].get()
		self.remove_popup(popup)
		
		if len( self.canvas.find_withtag('controls') ) == 4:
			self.redraw_handle_affected_curve()
			self.redraw_approximated_bezier()	
		return
		
	def delete_handle(self, id, popup):
	
		self.canvas.delete(id)
		self.remove_popup(popup)
		
		if len( self.canvas.find_withtag('controls') ) == 4:
			self.redraw_handle_affected_curve()
			self.redraw_approximated_bezier()	
		return	
		
	def remove_popup(self, popup):
	
		self.canvas.delete(popup)

def main():
  
	root = Tk()
#	root.geometry("600x400+20+20")
	app = Window(root)
	root.mainloop()	 

if __name__ == '__main__': main()
