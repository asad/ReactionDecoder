/*
 * Copyright (C) 2003-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package uk.ac.ebi.reactionblast.mapping.graph;

import static java.lang.Thread.currentThread;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class Lock {

    private boolean isLocked;
    private Thread lockedBy;
    private int lockedCount;

    Lock() {
        isLocked = false;
        lockedBy = null;
        lockedCount = 0;
    }

    public synchronized void lock() throws InterruptedException {
        Thread callingThread = currentThread();
        while (isIsLocked() && getLockedBy() != callingThread) {
            wait();
        }
        setIsLocked(true);
        setLockedCount(getLockedCount() + 1);
        setLockedBy(callingThread);
    }

    public synchronized void unlock() {
        if (currentThread() == this.getLockedBy()) {
            setLockedCount(getLockedCount() - 1);
            if (getLockedCount() == 0) {
                setIsLocked(false);
                notify();
            }
        }
    }

    /**
     * @return the isLocked
     */
    public boolean isIsLocked() {
        return isLocked;
    }

    /**
     * @param isLocked the isLocked to set
     */
    public void setIsLocked(boolean isLocked) {
        this.isLocked = isLocked;
    }

    /**
     * @return the lockedBy
     */
    public Thread getLockedBy() {
        return lockedBy;
    }

    /**
     * @param lockedBy the lockedBy to set
     */
    public void setLockedBy(Thread lockedBy) {
        this.lockedBy = lockedBy;
    }

    /**
     * @return the lockedCount
     */
    public int getLockedCount() {
        return lockedCount;
    }

    /**
     * @param lockedCount the lockedCount to set
     */
    public void setLockedCount(int lockedCount) {
        this.lockedCount = lockedCount;
    }
}
