#pragma once
//#include

#ifndef Q_MOC_RUN
#include "GLViewer.hpp"
#endif

#include <QObject>
#include <mutex>

namespace yade { // Cannot have #include directive inside.

/*
Singleton class managing OpenGL views,
a renderer instance and timer to refresh the display.
*/
class OpenGLManager : public QObject {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wsuggest-override"
	Q_OBJECT
#pragma GCC diagnostic pop
	DECLARE_LOGGER;

public:
	static OpenGLManager* self;
	OpenGLManager(QObject* parent = 0);
	// manipulation must lock viewsMutex!
	std::vector<shared_ptr<GLViewer>> views;
	std::vector<shared_ptr<GLViewer>> viewsLater;
	shared_ptr<OpenGLRenderer>        renderer;
	// signals are protected, emitting them is therefore wrapped with such funcs
	// NOTE: QT uses #define foreach Q_FOREACH which breaks boost. So I had to replace
	//       signals → Q_SIGNALS, slots → Q_SLOTS, emit → Q_EMIT also
	void emitResizeView(int id, int wd, int ht) { Q_EMIT resizeView(id, wd, ht); }
	void emitCreateView() { Q_EMIT createView(); }
	void emitStartTimer() { Q_EMIT startTimerSignal(); }
	void emitCloseView(int id) { Q_EMIT closeView(id); }
	// create a new view and wait for it to become available; return the view number
	// if timout (in seconds) elapses without the view to come up, reports error and returns -1
	int waitForNewView(double timeout = 5., bool center = true);
	// for commands yade.qt.center(…) and yade.qt.centerValues()
	Real     getSuggestedRadius() const;
	Vector3r getGridOrigin() const;
	Vector3r getSuggestedCenter() const;
	int      getGridDecimalPlaces() const;
Q_SIGNALS:
	void createView();
	void resizeView(int id, int wd, int ht);
	void closeView(int id);
	// this is used to start timer from the main thread via postEvent (ugly)
	void startTimerSignal();
public Q_SLOTS:
	virtual void cleanupViewsSlot();
	virtual void createViewSlot();
	virtual void resizeViewSlot(int id, int wd, int ht);
	virtual void closeViewSlot(int id = -1);
	void         timerEvent(QTimerEvent* event) override;
	virtual void startTimerSlot();
	void         centerAllViews(const Real& suggestedRadius, const Vector3r& gridOrigin, const Vector3r& suggestedCenter, int gridDecimalPlaces);

private:
	mutable std::mutex viewsMutex;
};

} // namespace yade
