/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for forming SCHP5 packets
 * @version 0.1
 * @date 2025-10-02
 */
# ifndef PACKET_FORMER_H
# define PACKET_FORMER_H

#include "packet.h"
#include "ring.hpp"

/**
 * @brief Forms packets according to the SCHP5 protocol for packet radios
 */
class PacketFormer {

    private:

    /**
     * @brief Contains iterator to current position in input buffer
     */
    std::ring<uint8_t>::iterator inputPosition;

    /**
     * @brief Contains iterator to current position in output buffer
     */
    std::ring<uint8_t>::iterator outputPosition;

    std::ring<uint8_t> test;

    public:

    /**
     * @brief Pointer to the buffer where the packet former reads in arbitrary binary data formatted as a buffer of bytes
     */
    std::ring<uint8_t> * inputBuffer;

    /**
     * @brief Pointer to the buffer where the formed packets are stored
     */
    std::ring<Packet> * outputBuffer;

    /**
     * @brief Construct a new packet former
     */
    PacketFormer(std::ring<uint8_t> * input_buffer, std::ring<Packet> * output_buffer);

};

#endif