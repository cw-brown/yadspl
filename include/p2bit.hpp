/**
 * @file p2bit.hpp
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines methods for converting SCHP5 packets to bits
 * @version 0.1
 * @date 2025-11-13
 */

#include "packet.h"

/**
 * @brief Writes a single packet to the data buffer
 * @param to_write the packet to write
 * @param position a pointer to the position to write the data
 * @returns the number of bits that have been written
 */
size_t singlePacketToBuffer(const Packet to_write, bool * position) {

    bool * start;

    start = position;

    // Write the first 8 flag bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.flag & (0b1 << (7 - i));
        position++;

    }

    // Write the 4 control code bits
    for(int i = 0; i < 4; i++) {

        *position = to_write.controlCode & (0b1 << (3 - i));
        position++;

    }

    // Write the 4 data length code bits
    for(int i = 0; i < 4; i++) {

        *position = to_write.dataLength & (0b1 << (3 - i));
        position++;

    }

    // Write the 8 receiver address bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.recieverAddress & (0b1 << (7 - i));
        position++;

    }

    // Write the 8 sender address bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.senderAddress & (0b1 << (7 - i));
        position++;

    }

    // Write the 8 sequence bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.sequenceNumber & (0b1 << (7 - i));
        position++;

    }

    // Write the data bits
    for(int i = 0; i < to_write.dataLength; i++) {

        for(int j = 0; j < 8; j++) {

            *position = to_write.data.at(i) & (0b1 << (7 - i));
            position++;

        }

    }

    // Write the 32 erc bits
    for(int i = 0; i < 32; i++) {

        *position = to_write.erc & (0b1 << (7 - i));
        position++;

    }

    // Write the final 8 flag bits
    for(int i = 0; i < 8; i++) {

        *position = to_write.flag & (0b1 << (7 - i));
        position++;

    }

    // Return total number of bits written
    return std::distance(start, position);

}